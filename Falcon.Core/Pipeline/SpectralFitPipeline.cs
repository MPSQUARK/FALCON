using Falcon.Core.Abstractions;
using Falcon.Core.Configuration;
using Falcon.Core.Diagnostics;
using Falcon.Core.Constants;
using Falcon.Core.Numerics;
using Falcon.Core.DependencyInjection;
using Falcon.Core.Services.Cosmology;
using Falcon.Core.Services.IO;
using Falcon.Core.Services.Fitting;
using Falcon.Core.Services.Instrument;
using Falcon.Core.Services.Statistics;
using ILGPU;
using ILGPU.Algorithms;
using ILGPU.Runtime;
using BAVCL;
using BAVCL.Core;
using BAVCL.Modules.Masking;
using BAVCL.Modules.Statistics;
using BAVCL.Types;
using System.Runtime.CompilerServices;

namespace Falcon.Core.Pipeline;

public class SpectralFitPipeline
{
    sealed class PipelineGpuKernels(
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<int>, int> trimDownFlux,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float> consecutive2D,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int> chiCalcRun1,
        Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>, ArrayView<int>, int> chiCalcRun2)
    {
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<int>, int> TrimDownFlux { get; } = trimDownFlux;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float> Consecutive2D { get; } = consecutive2D;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int> ChiCalcRun1 { get; } = chiCalcRun1;
        public Action<AcceleratorStream, Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>, ArrayView<int>, int> ChiCalcRun2 { get; } = chiCalcRun2;
    }

    static readonly ConditionalWeakTable<GPU, PipelineGpuKernels> PipelineKernelCache = new();

    static PipelineGpuKernels GetPipelineKernels(GPU gpu) =>
        PipelineKernelCache.GetValue(gpu, static g =>
        {
            var acc = g.accelerator;
            return new PipelineGpuKernels(
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<int>, int>(TrimDownFlux),
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, float, float>(ConsecutiveC_VDSP_2DKernel),
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int>(ChiCalcRun1),
                acc.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<int>, ArrayView<int>, int>(ChiCalcRun2));
        });

    readonly PipelineGpuKernels _pipelineKernels;
    public SpectralFitPipeline(GPU gpu,
                   string path,
                   bool milky_Way_Reddening = true,
                   bool hPF_Mode = true,
                   ushort n_Masked_Amstrongs = 20,
                   IDustService? dustService = null,
                   IStatisticsService? statisticsService = null,
                   ISpectralFitter? spectralFitter = null,
                   IFireflyFitsWriter? fitsWriter = null,
                   CancellationToken cancellationToken = default)
    {
        this.gpu = gpu ?? throw new ArgumentNullException(nameof(gpu));
        FalconGpuHost.EnsureReady(this.gpu);
        _pipelineKernels = GetPipelineKernels(this.gpu);
        DataFilePath = path;
        Milky_Way_Reddening = milky_Way_Reddening;
        HPF_Mode = hPF_Mode;
        N_Masked_Amstrongs = n_Masked_Amstrongs;
        _cancellationToken = cancellationToken;
        _dustService = dustService ?? new Services.Dust.DustService(gpu);
        _statisticsService = statisticsService ?? new StatisticsService();
        _spectralFitter = spectralFitter ?? new IterativeSpectralFitter(gpu);
        _fitsWriter = fitsWriter ?? new FireflyFitsWriter();
    }

    readonly GPU gpu;
    readonly IDustService _dustService;
    readonly IStatisticsService _statisticsService;
    readonly ISpectralFitter _spectralFitter;
    readonly IFireflyFitsWriter _fitsWriter;
    readonly CancellationToken _cancellationToken;

    // SPECTRA SPECIFIC VARIABLES
    #region
    // constructor variables
    public string DataFilePath { get; private set; }
    public bool Milky_Way_Reddening { get; private set; }
    public bool HPF_Mode { get; private set; }
    public ushort N_Masked_Amstrongs { get; private set; }

    // Initialisation variables
    /// <summary>
    /// DATA wavelength in OBSERVERS Frame of REFERENCE
    /// </summary>
    public Vector Wavelength { get; private set; }
    /// <summary>
    /// DATA Flux Parameter
    /// </summary>
    public Vector Flux { get; private set; }
    /// <summary>
    /// DATA ERROR PARAMETER
    /// </summary>
    public Vector Error { get; private set; }
    public float Redshift { get; private set; }
    public float[] RA_DEC { get; private set; }
    public float Velocity_Dispersion { get; private set; }
    public Vector Instrument_Resolution { get; private set; }

    public float Distance_Luminosity { get; private set; }
    public Vector Bad_Flags { get; private set; }
    public Vector Restframe_Wavelength { get; private set; }
    // !!errorNES!!
    public byte Trust_Flag { get; private set; }
    // !!errorNES!!
    public byte ObjID { get; private set; }
    // Bad Data ?? 
    // mask emission lines
    public float Ebv_MW { get; private set; }
    public float Ebv_Fit { get; private set; }
    public bool DataLiesOutsideModel = false;
    #endregion

    // MODEL SPECIFIC VARIABLES
    #region
    int Velocity_dispersion_r { get; set; }
    //int fit_per_iteration_cap { get; set; }

    // Model Main Values
    /// <summary>
    /// MODEL wavelengths
    /// </summary>
    public Vector Model_wavelength;
    /// <summary>
    /// MODEL FLUXES - List of MODEL FLUXES
    /// </summary>
    public Vector Model_flux; // Becomes NORMED after Method - NormaliseSpec
    /// <summary>
    /// MODEL AGES
    /// </summary>
    public float[] Model_ages;
    /// <summary>
    /// MODEL METALICITIES
    /// </summary>
    public float[] Model_metals;
    public float[] Mass_factor;
    public int[] MatchingWavelengthMapping;
    float[]? _modelFluxCpu;

    #endregion

    // Spectrum Initialisation
    public void InitialiseSpectraParameters(BAVCL.Vector Data, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
    {

        Wavelength = Data.GetColumnAsVectorX(0);
        Flux = Data.GetColumnAsVectorX(1);
        Error = Data.GetColumnAsVectorX(2);

        this.Redshift = Redshift;
        this.RA_DEC = RA_DEC;
        Velocity_Dispersion = Velocity_Disp;
        Trust_Flag = 1;
        ObjID = 0;

        if (Milky_Way_Reddening)
        {
            Ebv_MW = GetDustRADEC(this.RA_DEC, "ebv");
        }
        else
        {
            Ebv_MW = 0f;
        }

        Restframe_Wavelength = Wavelength * (1 / (1 + this.Redshift));

        ApplyGoodDataMask();

        Instrument_Resolution = Vector.Fill(gpu, instrument_resolution, Wavelength.Length);

        // CALCULATE LUMINOSITY DISTANCE in CM and convert to Mpc
        Distance_Luminosity = FitsSpectrumReader.Mpc2cm(LegacyConfig.Cosmology.GPU_Integration(gpu, this.Redshift, 1e-8f));


        // Initialise Model Section
        InitialiseSPModel();

    }

    public void InitialiseSpectraParameters(Vector Wavelength, Vector Flux, Vector Error, float Redshift, float[] RA_DEC, float Velocity_Disp, float instrument_resolution) // also include emission lines masking
    {

        // Set Long tasks to execute
        this.Wavelength = Wavelength;

        // Set Long tasks to execute
        this.Flux = Flux;

        // Set Long tasks to execute
        this.Error = Error;

        this.Redshift = Redshift;
        this.RA_DEC = RA_DEC;
        Velocity_Dispersion = Velocity_Disp;
        Trust_Flag = 1;
        ObjID = 0;


        if (Milky_Way_Reddening)
        {
            Ebv_MW = GetDustRADEC(this.RA_DEC, "ebv");
        }
        else
        {
            Ebv_MW = 0f;
        }

        Restframe_Wavelength = this.Wavelength * XMath.Rcp(1 + this.Redshift);

        ApplyGoodDataMask();

        Instrument_Resolution = Vector.Fill(gpu, instrument_resolution, this.Wavelength.Length);

        // CALCULATE LUMINOSITY DISTANCE in CM and convert to Mpc
        Distance_Luminosity = FitsSpectrumReader.Mpc2cm(LegacyConfig.Cosmology.GPU_Integration(gpu, this.Redshift, 1e-8f));


        // Initialise Model Section
        InitialiseSPModel();

    }



    // Model Initialisation
    private void InitialiseSPModel()
    {
        Velocity_dispersion_r = (int)(MathF.Round(Velocity_Dispersion / 5f) * 5f);
    }


    // SPECTRAL - SPECTRUM FUNCTIONS
    #region
    private void ApplyGoodDataMask()
    {
        Mask goodMask = BuildGoodDataMask(Flux, Error);
        if (!MaskHasBadLanes(goodMask))
            return;

        Vector filteredFlux = Flux.FilterX(goodMask);
        if (filteredFlux.Length == 0)
            throw new InvalidOperationException("No valid spectral pixels remain after quality masking.");

        Wavelength = Wavelength.FilterX(goodMask);
        Restframe_Wavelength = Restframe_Wavelength.FilterX(goodMask);
        Flux = filteredFlux;
        Error = Error.FilterX(goodMask);
    }

    /// <summary>
    /// Quality mask built entirely on the GPU — avoids pulling flux/error to CPU.
    /// Matches firefly: finite, strictly positive flux and error.
    /// </summary>
    private static Mask BuildGoodDataMask(Vector flux, Vector error)
    {
        const float upper = float.MaxValue;
        return (flux > 0f) & (flux < upper) & (error > 0f) & (error < upper);
    }

    /// <summary>Checks only the compact mask bitmap (not the full spectrum vectors).</summary>
    private static bool MaskHasBadLanes(Mask mask)
    {
        foreach (bool bit in mask.ToBoolArray())
        {
            if (!bit)
                return true;
        }

        return false;
    }

    private float GetDustRADEC(float[] RADEC, string dustmap, bool interpolate = true)
    {

        float[] l_b = Eq2Gal(
            FitsSpectrumReader.Degree2Radians(RADEC[0]),
            FitsSpectrumReader.Degree2Radians(RADEC[1])
            );

        return Get_SFD_dust(l_b[0], l_b[1]);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private float[] Eq2Gal(float ra, float dec)
    {

        float Cos_dec_Mult_Cos_GNEV0 = (float)(Math.Cos(dec) * Math.Cos(ra - FireflyConstants.Galactic_Northpole_Equatorial[0]));

        float b = (float)Math.Asin(
            Math.Sin(dec) * Math.Sin(FireflyConstants.Galactic_Northpole_Equatorial[1]) +
            Math.Cos(FireflyConstants.Galactic_Northpole_Equatorial[1]) * Cos_dec_Mult_Cos_GNEV0
            );

        float l = (float)Math.Atan2(
            Math.Sin(dec) * Math.Cos(FireflyConstants.Galactic_Northpole_Equatorial[1]) -
            Math.Sin(FireflyConstants.Galactic_Northpole_Equatorial[1]) * Cos_dec_Mult_Cos_GNEV0,
            Math.Cos(dec) * Math.Sin(ra - FireflyConstants.Galactic_Northpole_Equatorial[0])
            ) + FitsSpectrumReader.Degree2Radians(33f);

        if (l < 0)
        {
            l += FireflyConstants.TwoPi;
        }

        l %= FireflyConstants.TwoPi;


        return new float[2] { l * FireflyConstants.Rad2DegFactor, b * FireflyConstants.Rad2DegFactor };
    }

    private float Get_SFD_dust(float logitude, float latitude, string dustmap = "ebv", bool interpolate = true)
    {
        //Console.WriteLine($"longitude {logitude} , latitude {latitude}");

        // get path to dust map

        // if var dustmap is NOT a string 
        //      - Raise Exception
        //  + ALWAYS WILL BE A STRING

        // var dml = dustmap in lowercase

        /* if dustmap is 'ebv' 'eb-b' or 'e(b-v)'
         *  + ONLY accepted variation will be 'ebv'
         * 
         *      - read in file at os.environ['FF_DIR']+'/data/SFD_dust_4096_%s.fits'
         * + TEST WHERE THIS IS TO CONVERT TO HDF5
         * 
         * + RETURN ALL OTHER TYPES AS DUSTMAP NOT FOUND
         * 
         * if dustmap = i100
         * if dustmap = x
         * if dustmap = t
         * if dustmap = mask
         * else dustmapfn = dustmap ?? would that not crash?
         * 
         * if longitude is a scalar
         *      - l = convert longitude to an array and multipy by pi/180
         * + WHICH it IS
         * if latitude is a scalar
         *      - b = convert latitude to an array and multipy by pi/180
         *      
         * if length of l and b DONT match, RAISE EXCEPTION
         *  + WILL ALWAYS MATCH SINCE l and b are floats
         * 
         * if NOT dustmapfn.contains(%s)
         *  + %s WILL ALWAYS be in the path name
         *      - f = fits.open(dustmapfn)
         *      - try : mapds = [f[0].data] then f.close()
         * + READ IN DUSTMAP in python and convert to HDF5 for this and read in PreInitialise Phase
         * 
         *      - assert that mapds[-1].shape[0] is = to mapds[-1].shape[1] 
         * + Test in Python what does means
         * 
         *      - polename = dustmapfn.split('.')[0].split('_')[-1].lower() ?? Can do a Contains or Regex test instead??
         * 
         *      - if polename is 'ngp'
         *          -> n=[1]
         *      - if sum( b > 0 ) > 0
         *          -> b = b
         *      - if polename is 'sgp'
         *          -> n = [-1]
         *          - if sum ( b < 0 ) > 0:
         *              -> b = b
         *      - else
         *          - RAISE EXCEPTION
         *      
         *      -> masks = [ones(same shape as b) converted to boolean ]
         * ELSE
         *      -> nmask = b >= 0 
         *      -> smask = !nmask
         *      
         *      masks = [nmask, smask]
         *      ns = [1,-1]
         *      
         *      mapds = []
         *     =>-> READ IN THE dustmapfn fits for 'ngp'/'sgp'
         *      try: mapds.add(f[0].data) then close(f)
         *      asstert that mapds[-1].shape[0] is EQUAL to mapds[-1].shape[1]
         *     <=->repeat for 'sgp'
         *      
         *      retvals = []
         *      for n,mapd,m in zip(ns,mapds,masks) 
         *          + CHECK THIS FUNCTIONALITY IN PYTHON TO REPLICATE
         *          npix = mapd.shape[0]
         *      
         *      
         *      
        */



        return (float)_dustService.GetSfdEbv(logitude, latitude, FireflyConstants.ngp_dust, FireflyConstants.sgp_dust);
    }

    #endregion

    // MODELS - MODEL FUNCTIONS
    #region
    public void Fit_models_to_data()
    {
        Get_Model();

        if (LegacyConfig.FittingEngine == "legacy_gpu")
        {
            RunLegacyGpuPath();
            return;
        }

        var orchestrator = new FireflyFitOrchestrator(_dustService, _statisticsService, _spectralFitter, _fitsWriter);
        var ctx = BuildFireflyContext();
        if (orchestrator.TryFit(ctx, out _))
            return;

        if (LegacyConfig.LegacyGpuFallback)
            RunLegacyGpuPath();
    }

    private void RunLegacyGpuPath()
    {
        NormaliseSpec();
        using (GpuScope.BeginReadOnly(Model_flux, Flux, Error))
        {
            if (LegacyConfig.HPF_Mode)
                RunDustAttenuationPipeline();
            GenerateInterpSSP();
        }
    }

    private static string DustLawName() => LegacyConfig.Dust_Law switch
    {
        1 => "allen",
        2 => "prevot",
        _ => "calzetti"
    };

    private FireflyFitContext BuildFireflyContext()
    {
        int nModels = Model_ages.Length;
        int nData = Restframe_Wavelength.Length;
        int cols = Model_flux.Columns;
        double[] flat = _modelFluxCpu is { Length: > 0 } cpu
            ? Array.ConvertAll(cpu, v => (double)v)
            : Model_flux.RetrieveReadOnlySpan().ToArray().Select(v => (double)v).ToArray();
        int[] map = MatchingWavelengthMapping ?? Enumerable.Range(0, Math.Min(nData, cols)).ToArray();

        var modelFluxRaw = new double[nModels][];
        for (int m = 0; m < nModels; m++)
        {
            modelFluxRaw[m] = new double[nData];
            int rowOffset = m * cols;
            for (int i = 0; i < nData; i++)
            {
                int modelIdx = i < map.Length ? map[i] : Math.Min(i, cols - 1);
                modelIdx = Math.Clamp(modelIdx, 0, cols - 1);
                double v = flat[rowOffset + modelIdx];
                modelFluxRaw[m][i] = double.IsFinite(v) ? v : 0;
            }
        }

        ReadOnlySpan<float> dataWave = Restframe_Wavelength.RetrieveReadOnlySpan();
        ReadOnlySpan<float> dataFlux = Flux.RetrieveReadOnlySpan();
        ReadOnlySpan<float> dataError = Error.RetrieveReadOnlySpan();

        var dataWaveArr = new double[nData];
        var dataFluxArr = new double[nData];
        var dataErrorArr = new double[nData];
        for (int i = 0; i < nData; i++)
        {
            dataWaveArr[i] = dataWave[i];
            dataFluxArr[i] = dataFlux[i];
            dataErrorArr[i] = dataError[i];
        }

        return new FireflyFitContext
        {
            DataFilePath = DataFilePath,
            DataWavelength = dataWaveArr,
            DataFlux = dataFluxArr,
            DataError = dataErrorArr,
            ModelWavelength = dataWaveArr,
            ModelFluxRaw = modelFluxRaw,
            ModelsOnDataGrid = true,
            ModelAges = Model_ages.Select(a => (double)a).ToArray(),
            ModelMetals = Model_metals.Select(z => (double)z).ToArray(),
            DistanceLumCm = Distance_Luminosity,
            Redshift = Redshift
        };
    }

    static void TraceGetModelComplete(int nModels, float[] fluxBank, float[] ages)
    {
        if (!FitTrace.Enabled) return;
        int cols = nModels > 0 ? fluxBank.Length / nModels : 0;
        var model0 = cols > 0 ? fluxBank.AsSpan(0, cols).ToArray() : [];
        double bankMax = 0;
        int bankPositive = 0;
        foreach (double v in fluxBank)
        {
            if (double.IsFinite(v) && v > 0) bankPositive++;
            if (double.IsFinite(v) && v > bankMax) bankMax = v;
        }
        FitTrace.Emit("get_model", new Dictionary<string, object>
        {
            ["n_models"] = nModels,
            ["n_wave"] = cols,
            ["age_min"] = ages.Length > 0 ? ages.Min() : 0,
            ["age_max"] = ages.Length > 0 ? ages.Max() : 0,
            ["bank_max"] = bankMax,
            ["bank_positive"] = bankPositive,
        }.Concat(FitTrace.Stats("raw_model0", model0.Select(v => (double)v).ToArray()))
         .ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    private void Get_Model()
    {
        using var _ = FitTrace.BeginSpan("get_model");
        // M11 Models
        if (LegacyConfig.Model_Key % 2 == 1)
        {
            //var first_file = true;
            // List of model_files = []

            // Constructs the metallicity array of models:
            //var all_metal_files = sorted(glob.glob(model_path+"*"))
            // metal_files = []
            // metal = [] #[-2.25,-1.35,-0.33,0,0.35]
            // for z in range(len(all_metal_files)):
            //.
            //.
            //.

            return;
        }

        // MaStar Models
        if (LegacyConfig.Model_Key % 2 == 0)
        {
            // Gets Wavelengths from Model Data
            Model_wavelength = FireflyConstants.wavelength;
            // Gets Indexes of Models matching wavlengths
            TrimModel();

            int nWave = FireflyConstants.wavelengthCpu.Length;
            double[]? mwAtten = null;
            if (Ebv_MW != 0)
            {
                var waveDoubles = new double[nWave];
                for (int k = 0; k < nWave; k++)
                    waveDoubles[k] = FireflyConstants.wavelengthCpu[k];
                mwAtten = _dustService.FitzpatrickUnred(waveDoubles, -Ebv_MW);
                if (FitTrace.Enabled)
                {
                    double aMin = double.PositiveInfinity, aMax = 0;
                    int aZero = 0;
                    foreach (double a in mwAtten)
                    {
                        if (a == 0) aZero++;
                        if (a < aMin) aMin = a;
                        if (a > aMax) aMax = a;
                    }
                    FitTrace.Emit("mw_atten", new Dictionary<string, object>
                    {
                        ["ebv_mw"] = Ebv_MW,
                        ["min"] = aMin,
                        ["max"] = aMax,
                        ["zeros"] = aZero,
                    });
                }
            }

            Vector modelSres = FireflyConstants.sres.Copy(Cache: false);
            (int[] indx, ResolutionService res, VariableGaussianKernel vGK) = DownGradeModelInvar(Model_wavelength, Velocity_dispersion_r, modelSres, Instrument_Resolution);

            List<double> model_flux = new();
            List<float> age_model = new();
            List<float> metal_model = new();
            List<(float[] Row, float Age, float Metal)> pendingRows = LegacyConfig.Downgrade_models ? new() : null!;

            float[] rowF = new float[nWave];

            for (int i = 0; i < FireflyConstants.t.Length; i++)
            {
                _cancellationToken.ThrowIfCancellationRequested();
                //index i in t , t value
                if ((FireflyConstants.t[i] < LegacyConfig.MinMax_Model_Age[0]) || (FireflyConstants.t[i] > LegacyConfig.MinMax_Model_Age[1]))
                {
                    continue;
                }

                for (int j = 0; j < FireflyConstants.Z.Length; j++)
                {
                    // index, Z value
                    if (FireflyConstants.Z[j] < LegacyConfig.MinMax_Model_Metalicity[0] || FireflyConstants.Z[j] > LegacyConfig.MinMax_Model_Metalicity[1])
                    {
                        continue;
                    }
                    if (FireflyConstants.Z[j] < -1.35f && FireflyConstants.t[i] < 1)
                    {
                        continue;
                    }

                    for (int k = 0; k < nWave; k++)
                        rowF[k] = SanitizeModelFlux(FireflyConstants.fluxgrid[i, j, FireflyConstants.sidx, k]);

                    if (LegacyConfig.Downgrade_models)
                    {
                        pendingRows.Add(((float[])rowF.Clone(), FireflyConstants.t[i], FireflyConstants.Z[j]));
                        continue;
                    }

                    AppendModelRow(model_flux, rowF, mwAtten, null, age_model, metal_model, FireflyConstants.t[i], FireflyConstants.Z[j]);
                }
            }

            if (LegacyConfig.Downgrade_models && pendingRows.Count > 0)
            {
                BatchMatchSpectralResolution(pendingRows, indx, vGK);
                var mwScratch = mwAtten is null ? null : new double[nWave];
                foreach (var (row, age, metal) in pendingRows)
                    AppendModelRow(model_flux, row, mwAtten, mwScratch, age_model, metal_model, age, metal);
            }

            _modelFluxCpu = model_flux.Select(v => v > float.MaxValue ? float.MaxValue : (float)v).ToArray();
            var gpuFlux = (float[])_modelFluxCpu.Clone();
            Model_flux = new Vector(gpu, gpuFlux, nWave);
            Model_ages = age_model.ToArray();
            Model_metals = metal_model.ToArray();
            TraceGetModelComplete(Model_ages.Length, _modelFluxCpu, Model_ages);

            return;
        }


    }

    private void TrimModel()
    {
        ReadOnlySpan<float> modelWave = Model_wavelength.RetrieveReadOnlySpan();
        ReadOnlySpan<float> dataWave = Restframe_Wavelength.RetrieveReadOnlySpan();

        int length_Data = dataWave.Length;
        int length_Mod = modelWave.Length;
        int[] indxs;

        if (modelWave[0] < dataWave[0] && modelWave[^1] > dataWave[^1])
        {
            int endIdx = IndexOfMinAbsDiff(modelWave, dataWave[^1]) + 1;
            int startIdx = IndexOfMinAbsDiff(modelWave[..(endIdx + 1)], dataWave[0]);

            if (endIdx - startIdx == length_Data)
            {
                MatchingWavelengthMapping = Enumerable.Range(startIdx, endIdx - startIdx).ToArray();
                return;
            }

            if (endIdx - startIdx != length_Data)
            {
                indxs = new int[length_Data];
                indxs[0] = startIdx;
                indxs[^1] = endIdx;
                for (int i = 1; i < indxs.Length - 1; i++)
                    indxs[i] = IndexOfMinAbsDiff(modelWave[..endIdx], dataWave[i]);

                MatchingWavelengthMapping = indxs;
                return;
            }
        }

        DataLiesOutsideModel = true;

        indxs = new int[length_Data];
        for (int i = 0; i < dataWave.Length; i++)
            indxs[i] = IndexOfNearestSorted(modelWave, dataWave[i]);

        MatchingWavelengthMapping = indxs;
    }

    static int IndexOfNearestSorted(ReadOnlySpan<float> sortedValues, float target)
    {
        if (sortedValues.IsEmpty) return 0;
        int lo = 0;
        int hi = sortedValues.Length - 1;
        while (lo < hi)
        {
            int mid = (lo + hi) >> 1;
            if (sortedValues[mid] < target) lo = mid + 1;
            else hi = mid;
        }
        if (lo > 0 && MathF.Abs(sortedValues[lo - 1] - target) < MathF.Abs(sortedValues[lo] - target))
            return lo - 1;
        return lo;
    }

    /// <summary>MaStar fluxgrid uses -99 as missing-pixel sentinel (firefly leaves other values unchanged).</summary>
    private static float SanitizeModelFlux(float value)
    {
        if (!float.IsFinite(value) || value <= -90f)
            return 0f;
        return value;
    }

    private static int IndexOfMinAbsDiff(ReadOnlySpan<float> values, float target)
    {
        int best = 0;
        float bestDiff = MathF.Abs(values[0] - target);
        for (int i = 1; i < values.Length; i++)
        {
            float diff = MathF.Abs(values[i] - target);
            if (diff < bestDiff)
            {
                bestDiff = diff;
                best = i;
            }
        }

        return best;
    }

    private static Vector MultiplyElementWiseMatchingLength(Vector left, Vector right)
    {
        if (left.Shape().MatchesDimensions(right.Shape()))
            return left * right;

        int leftCols = left.Columns;
        int rightCols = right.Columns;
        left.Flatten();
        right.Flatten();
        Vector product = left * right;
        left.Columns = leftCols;
        right.Columns = rightCols;
        return product;
    }

    private (int[], ResolutionService, VariableGaussianKernel) DownGradeModelInvar(Vector mod_wavelength, int vdisp_round, Vector sres, Vector r_instrument, float min_sig_pix = 0f, bool no_offset = true)
    {
        Vector new_sig = new(gpu, new float[mod_wavelength.Length], 0, false);
        ReadOnlySpan<float> rInst = r_instrument.RetrieveReadOnlySpan();
        var mappingHead = MatchingWavelengthMapping.Length > 1
            ? MatchingWavelengthMapping[..^1].ToHashSet()
            : new HashSet<int>();

        using (var scope = new_sig.CpuScopeAndSync())
        {
            EditableView<float> view = scope.View;
            for (int i = 0, j = 0; i < mod_wavelength.Length; i++)
            {
                if (mappingHead.Contains(i))
                    j++;

                float sig_instrument = FireflyConstants.c_kms / (rInst[j] * FireflyConstants.sig2FWHM);
                view[i] = MathF.Sqrt(MathF.Pow(vdisp_round, 2f) + MathF.Pow(sig_instrument, 2f));
            }
        }

        Vector new_sres = FireflyConstants.c_div_sig2 / new_sig;

        bool log_wave = true;
        float step12 = mod_wavelength.GetAt(3) - mod_wavelength.GetAt(2);
        float step01 = mod_wavelength.GetAt(2) - mod_wavelength.GetAt(1);
        if (step12 - step01 < 0.000001f * step01)
            log_wave = false;


        // match_ResolutionService Begins

        Vector wave = mod_wavelength;
        Vector new_sres_wave = mod_wavelength;

        if (wave.MinX() < new_sres_wave.GetAt(0) || wave.MaxX() > new_sres_wave.GetAt(new_sres_wave.Length - 1))
        {
            Console.WriteLine("WARNING : Mapping to the new spectral resolution will require extrapolating the provided input vectors!");
        }

        ResolutionService new_res = new(gpu, new_sres_wave, new_sres, log_wave);
        ResolutionService res = new(gpu, wave, sres, log_wave);

        res.match(new_res, no_offset, min_sig_pix);

        Mask activeMask = res.sig_pd > min_sig_pix;
        bool[] activeBits = activeMask.ToBoolArray();
        var indxList = new List<int>(activeBits.Length);
        for (int i = 0; i < activeBits.Length; i++)
        {
            if (activeBits[i])
                indxList.Add(i);
        }

        ReadOnlySpan<float> sigPd = res.sig_pd.RetrieveReadOnlySpan();
        int[] indx = indxList.ToArray();
        int len = indx.Length;
        VariableGaussianKernel vGK;

        if (len > 0f)
        {
            float[] selectedsig = new float[len];
            for (int i = 0; i < len; i++)
                selectedsig[i] = sigPd[indx[i]];

            vGK = new VariableGaussianKernel(gpu, selectedsig);
        }
        else
        {
            vGK = new VariableGaussianKernel(gpu, res.sig_pd.RetrieveReadOnlySpan().ToArray());
        }

        //if (res.sig_mask.Contains(true) || mask.Contains(true))
        //{
        //    for (int i = 0; i < mask.Length; i++)
        //    {
        //        mask[i] = res.sig_mask[i] || mask[i];  // mask_out
        //    }
        //}

        return (indx, res, vGK);
    }

    static void AppendModelRow(
        List<double> modelFlux,
        ReadOnlySpan<float> rowF,
        double[]? mwAtten,
        double[]? mwScratch,
        List<float> ageModel,
        List<float> metalModel,
        float age,
        float logMetal)
    {
        if (mwAtten is null)
        {
            for (int k = 0; k < rowF.Length; k++)
            {
                float f = SimdSpanOps.SanitizeFluxScalar(rowF[k]);
                modelFlux.Add(f);
            }
        }
        else
        {
            if (mwScratch is null || mwScratch.Length < rowF.Length)
                throw new ArgumentException("mwScratch must cover row length when Milky Way attenuation is enabled.", nameof(mwScratch));
            SimdSpanOps.MultiplyFloatToDouble(rowF, mwAtten, mwScratch.AsSpan(0, rowF.Length));
            for (int k = 0; k < rowF.Length; k++)
                modelFlux.Add(mwScratch[k]);
        }

        ageModel.Add(age);
        metalModel.Add(MathF.Pow(10f, logMetal));
    }

    void BatchMatchSpectralResolution(List<(float[] Row, float Age, float Metal)> rows, int[] indx, VariableGaussianKernel vGK)
    {
        for (int r = 0; r < rows.Count; r++)
            MatchSpectralResolution(gpu, rows[r].Row, indx, vGK);
    }

    private float[] DownGrade(GPU gpu, float[] flux, int[] indx, ResolutionService res, VariableGaussianKernel vGK)
    {
        MatchSpectralResolution(gpu, flux, indx, vGK);
        return flux;
    }

    private static void MatchSpectralResolution(GPU gpu, float[] flux, int[] indx, VariableGaussianKernel vGK)
    {
        float peak = 0f;
        for (int i = 0; i < flux.Length; i++)
            if (float.IsFinite(flux[i]) && flux[i] > peak) peak = flux[i];
        if (peak <= 0f) return;

        float invPeak = 1f / peak;
        SimdSpanOps.ScaleInPlace(flux, invPeak);

        int len = indx.Length;
        if (len > 0)
        {
            float[] selectedflux = new float[len];
            for (int i = 0; i < len; i++)
                selectedflux[i] = flux[indx[i]];

            float[] selected_outflux = vGK.Convolve(gpu, selectedflux);
            for (int i = 0; i < len; i++)
                flux[indx[i]] = selected_outflux[i];
        }
        else
        {
            Array.Copy(vGK.Convolve(gpu, flux), flux, flux.Length);
        }

        SimdSpanOps.ScaleInPlace(flux, peak);
    }



    private void Match_data_models()
    {
        // Get model_flux and Data_flux




        /*
         * sets the
         * wave
         * model_flux_raw
         * 
         * outputs the 
         * data_flux
         * error_flux
         */
    }


    private void NormaliseSpec()
    {
        float data_norm;
        using (GpuScope.BeginReadOnly(Flux))
            data_norm = Flux.MedianX();

        int num_mods = Model_flux.Length / Model_flux.Columns;

        Vector model_norm = new(gpu, new float[num_mods], 0, false);
        float[] mass_factor = new float[num_mods];
        ReadOnlySpan<float> modelFlux;
        using (GpuScope.BeginReadOnly(Model_flux))
            modelFlux = Model_flux.RetrieveReadOnlySpan();
        int cols = Model_flux.Columns;

        using (var scope = model_norm.CpuScopeAndSync())
        {
            EditableView<float> view = scope.View;
            for (int i = 0; i < num_mods; i++)
                view[i] = MedianOfSpan(modelFlux.Slice(i * cols, cols));
        }

        using (GpuScope.BeginReadOnly(Model_flux))
            Model_flux = ConsecutiveCompoundOperation2D(Model_flux, model_norm, data_norm);
        for (int i = 0; i < num_mods; i++)
            mass_factor[i] = data_norm / model_norm.RetrieveReadOnlySpan()[i];
        Mass_factor = mass_factor;
    }

    static float MedianOfSpan(ReadOnlySpan<float> values)
    {
        var scratch = new float[values.Length];
        values.CopyTo(scratch);
        return FitsSpectrumReader.Median(scratch);
    }

    public Vector ConsecutiveCompoundOperation2D(Vector vectorA, Vector vectorB, float scalar)
    {
        Vector output = new(gpu, new float[vectorA.Length], vectorA.Columns);

        using (GpuScope.Begin(output, vectorA, vectorB))
        {
            var buffer = output.GetBuffer();
            var buffer2 = vectorA.GetBuffer();
            var buffer3 = vectorB.GetBuffer();
            _pipelineKernels.Consecutive2D(
                gpu.accelerator.DefaultStream,
                buffer.IntExtent,
                buffer.View,
                buffer2.View,
                buffer3.View,
                vectorA.Columns,
                scalar);
            gpu.accelerator.Synchronize();
        }

        return output;
    }
    // KERNEL
    static void ConsecutiveC_VDSP_2DKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputA, ArrayView<float> InputB, float Cols, float Scalar)
    {
        Output[index] = (InputA[index] / InputB[(int)(index / Cols)]) * Scalar;
    }

    private Vector unred(Vector wavelength, float ebv_mw)
    {
        var wave = wavelength.RetrieveReadOnlySpan().ToArray();
        var atten = _dustService.FitzpatrickUnred(wave.Select(w => (double)w).ToArray(), ebv_mw);
        return new Vector(gpu, atten.Select(a => (float)a).ToArray());
    }

    private void RunDustAttenuationPipeline()
    {
        int nModels = Model_ages.Length;
        int nWave = Flux.Length;
        var wave = Restframe_Wavelength.RetrieveReadOnlySpan().ToArray().Select(w => (double)w).ToArray();
        var data = Flux.RetrieveReadOnlySpan().ToArray().Select(v => (double)v).ToArray();
        var err = Error.RetrieveReadOnlySpan().ToArray().Select(v => (double)v).ToArray();
        var modelFlux = new double[nModels][];
        for (int m = 0; m < nModels; m++)
            modelFlux[m] = Model_flux.GetRowAsArray(m).Select(v => (double)v).ToArray();

        var dustInput = new DustAttenuationInput
        {
            Wavelength = wave,
            DataFlux = data,
            ErrorFlux = err,
            ModelFlux = modelFlux,
            Age = Model_ages.Select(a => (double)a).ToArray(),
            Metal = Model_metals.Select(z => (double)z).ToArray(),
            MaxEbv = 1f,
            NumDustVals = 200,
            SmoothingLength = 100f,
            DustLaw = LegacyConfig.Dust_Law switch { 1 => "allen", 2 => "prevot", _ => "calzetti" }
        };

        var fitter = _spectralFitter;
        var result = _dustService.DetermineAttenuation(dustInput, fitter);
        Ebv_Fit = (float)result.Ebv;
    }


    #endregion


    // Test CODE
    public float CalculateChiSqu(int model)
    {
        int length = Restframe_Wavelength.Length;

        Vector absDiff = (Restframe_Wavelength.GetAt(0) - Model_wavelength).AbsXIP();
        int idx_closest = IndexOfMin(absDiff.RetrieveReadOnlySpan());

        ReadOnlySpan<float> modelFluxRow = Model_flux.RetrieveReadOnlySpan()
            .Slice((model * Model_flux.Columns) + idx_closest, length);
        ReadOnlySpan<float> flux = Flux.RetrieveReadOnlySpan();
        ReadOnlySpan<float> error = Error.RetrieveReadOnlySpan();

        float sum = 0f;
        for (int j = 0; j < length; j++)
            sum += MathF.Pow((modelFluxRow[j] - flux[j]) / error[j], 2f);

        return sum;
    }

    // will combine together SSP 50:50, 30:30:30, 25:25:25:25, 20:20:20:20:20
    private void GenerateInterpSSP()
    {
        int totlen = Model_flux.Length;
        int lenOneModel = Model_flux.Columns;
        int noModels = totlen / lenOneModel;

        using (GpuScope.BeginReadOnly(Model_flux, Flux, Error))
        {
            var trimKernel = _pipelineKernels.TrimDownFlux;
            Vector trimmedflux = new(gpu, new float[MatchingWavelengthMapping.Length * noModels], MatchingWavelengthMapping.Length);

            using (GpuScope.Begin(trimmedflux, Model_flux))
            using (MemoryBuffer1D<int, Stride1D.Dense> mappingBuffer = gpu.accelerator.Allocate1D(MatchingWavelengthMapping))
            {
                trimKernel(
                    gpu.accelerator.DefaultStream,
                    noModels,
                    trimmedflux.GetBuffer().View,
                    Model_flux.GetBuffer().View,
                    mappingBuffer.View,
                    lenOneModel);
                gpu.accelerator.Synchronize();
            }

            var chiKernel1 = _pipelineKernels.ChiCalcRun1;
            Vector Chis = new(gpu, new float[noModels]);

            using (GpuScope.Begin(Chis, trimmedflux, Flux, Error))
            {
                chiKernel1(
                    gpu.accelerator.DefaultStream,
                    noModels,
                    Chis.GetBuffer().View,
                    Flux.GetBuffer().View,
                    Error.GetBuffer().View,
                    trimmedflux.GetBuffer().View,
                    Flux.Length);
                gpu.accelerator.Synchronize();
            }

            List<int> pairX = [], pairY = [];

            for (int j = noModels - 1, i = 0; j >= 1; j--, i++)
            {
                pairX.AddRange(Enumerable.Repeat(i, j));
                pairY.AddRange(Enumerable.Range(i + 1, j));
            }

            var chiKernel2 = _pipelineKernels.ChiCalcRun2;
            Vector ChisG1 = new(gpu, new float[pairX.Count]);

            using (GpuScope.Begin(ChisG1, trimmedflux, Flux, Error))
            using (MemoryBuffer1D<int, Stride1D.Dense> bufferX = gpu.accelerator.Allocate1D(pairX.ToArray()))
            using (MemoryBuffer1D<int, Stride1D.Dense> bufferY = gpu.accelerator.Allocate1D(pairY.ToArray()))
            {
                chiKernel2(
                    gpu.accelerator.DefaultStream,
                    pairX.Count,
                    ChisG1.GetBuffer().View,
                    Flux.GetBuffer().View,
                    Error.GetBuffer().View,
                    trimmedflux.GetBuffer().View,
                    bufferX.View,
                    bufferY.View,
                    Flux.Length);
                gpu.accelerator.Synchronize();
            }

            float minChi = Chis.MinX();
            float minChiG1 = ChisG1.MinX();
            ReadOnlySpan<float> chisG1Cpu = ChisG1.RetrieveReadOnlySpan();
            int best_model = IndexOfMin(chisG1Cpu);
            int bestSingle = IndexOfMin(Chis.RetrieveReadOnlySpan());

            RecordData(bestSingle, minChi / Flux.Length, [pairX[best_model], pairY[best_model]], minChiG1 / Flux.Length);
        }
    }

    static int IndexOfMin(ReadOnlySpan<float> values)
    {
        int best = 0;
        float min = values[0];
        for (int i = 1; i < values.Length; i++)
        {
            if (values[i] < min)
            {
                min = values[i];
                best = i;
            }
        }

        return best;
    }

    // Custom ILGPU kernels: 2D model-bank χ² with per-branch clipped pixel lists (ChiCalcRun2).
    // BAVCL lacks an equivalent fused gather+reduce over variable index pairs per model row.
    static void TrimDownFlux(Index1D index, ArrayView<float> Output, ArrayView<float> Fluxes, ArrayView<int> idxs, int fluxlen)
    {
        for (int i = 0; i < idxs.Length; i++)
        {
            Output[index * idxs.Length + i] = Fluxes[index * fluxlen + idxs[i]];
        }
    }
    static void ChiCalcRun1(Index1D index, ArrayView<float> Chis, ArrayView<float> DataFlux, ArrayView<float> DataError, ArrayView<float> ModelFlux, int Modellen)
    {
        float sum = 0f;
        for (int i = 0; i < Modellen; i++)
        {
            sum += XMath.Pow((ModelFlux[index * Modellen + i] - DataFlux[i]) / DataError[i], 2f);
        }

        Chis[index] = sum;

    }
    static void ChiCalcRun2(Index1D index, ArrayView<float> Chis, ArrayView<float> DataFlux, ArrayView<float> DataError, ArrayView<float> ModelFlux, ArrayView<int> X, ArrayView<int> Y, int Modellen)
    {
        float sum = 0f;
        for (int i = 0; i < Modellen; i++)
        {
            sum += XMath.Pow((((ModelFlux[Modellen * X[index] + i] + ModelFlux[Modellen * Y[index] + i]) * 0.5f) - DataFlux[i]) / DataError[i], 2f);
        }

        Chis[index] = sum;
    }


    // record the data output
    private void RecordData(int SSPmodelnum, float SSPchi, int[] dualsspmodelnum, float dualsspchi)
    {
        string outDirectory = Path.Join(LegacyConfig.PathOfProgram, "/Output");
        if (!Directory.Exists(outDirectory))
            Directory.CreateDirectory(outDirectory);

        string data = DataFilePath.Replace("/Data", "/Output").Replace(".fits", "") + ".txt";
        File.WriteAllText(data, "");

        File.AppendAllText(data, "Restframe Wavelength of Data :\n");
        File.AppendAllText(data, string.Join(" ", Restframe_Wavelength.RetrieveReadOnlySpan().ToArray()));

        File.AppendAllText(data, "\nFlux of Data : \n");
        File.AppendAllText(data, string.Join(" ", Flux.RetrieveReadOnlySpan().ToArray()));

        File.AppendAllText(data, "\nError of Data : \n");
        File.AppendAllText(data, string.Join(" ", Error.RetrieveReadOnlySpan().ToArray()));

        File.AppendAllText(data, "\nModel Wavelengths : \n");
        File.AppendAllText(data, string.Join(" ", Model_wavelength.RetrieveReadOnlySpan().ToArray()));


        File.AppendAllText(data, "\n\nBest Single Model Fit : \n");
        File.AppendAllText(data, $"\nModel no : {SSPmodelnum} \n");
        File.AppendAllText(data, $"\nModel chi : {SSPchi} \n");
        File.AppendAllText(data, $"\nGalaxy Age Predicted : {Model_ages[SSPmodelnum]} Gyr \n");
        File.AppendAllText(data, $"\nGalaxy Metalicity Predicted : {Model_metals[SSPmodelnum]} [Z/H] \n");
        File.AppendAllText(data, $"\nGalaxy Mass Predicted : {(1f / Mass_factor[SSPmodelnum]) * 1e-17f} Msolar \n");

        File.AppendAllText(data, "\nModel Flux : \n");
        File.AppendAllText(data, string.Join(" ", Model_flux.GetRowAsArray(SSPmodelnum)));




        File.AppendAllText(data, "\n\nBest Dual-SSP Model Fit : \n");
        File.AppendAllText(data, $"\nModel no : {dualsspmodelnum[0]} & {dualsspmodelnum[1]} \n");
        File.AppendAllText(data, $"\nModel chi : {dualsspchi} \n");
        File.AppendAllText(data, $"\nGalaxy Age Predicted : 50% - {Model_ages[dualsspmodelnum[0]]} / 50% - {Model_ages[dualsspmodelnum[1]]} Gyr \n");
        File.AppendAllText(data, $"\nGalaxy Metalicity Predicted : {0.5f * (Model_metals[dualsspmodelnum[0]] + Model_metals[dualsspmodelnum[1]])} [Z/H] \n");
        File.AppendAllText(data, $"\nGalaxy Mass Predicted : {0.5f * (((1f / Mass_factor[dualsspmodelnum[0]]) * 1e-17f) + ((1f / Mass_factor[dualsspmodelnum[1]]) * 1e-17f))} Msolar \n");

        float[] midptflux = Midpoint(
            Model_flux.GetRowAsArray(dualsspmodelnum[0]),
            Model_flux.GetRowAsArray(dualsspmodelnum[1]));

        File.AppendAllText(data, "\nModel Flux : \n");
        File.AppendAllText(data, string.Join(" ", midptflux));

    }

    private static float[] Midpoint(float[] vectA, float[] vectB)
    {
        float[] vectR = new float[vectA.Length];
        for (int i = 0; i < vectA.Length; i++)
        {
            vectR[i] = 0.5f * (vectA[i] + vectB[i]);
        }
        return vectR;
    }



    // EXPERIMENTAL CODE DONT USE YET
    public float[] CalculateChiSquVec()
    {
        ReadOnlySpan<float> modelWave = Model_wavelength.RetrieveReadOnlySpan();
        ReadOnlySpan<float> dataWave = Restframe_Wavelength.RetrieveReadOnlySpan();
        int length_Data = dataWave.Length;
        int length_Mod = modelWave.Length;
        int models = length_Data / length_Mod;

        int endIdx = IndexOfMinAbsDiff(modelWave, dataWave[^1]) + 1;
        int startIdx = IndexOfMinAbsDiff(modelWave[..(endIdx + 1)], dataWave[0]);
        int model_len = endIdx - startIdx;

        ReadOnlySpan<float> modelFlux = Model_flux.RetrieveReadOnlySpan();
        List<float> new_mod_flux = new(models * model_len);
        for (int i = 0; i < models; i++)
        {
            int rowStart = i * length_Mod + startIdx;
            for (int j = 0; j < model_len; j++)
                new_mod_flux.Add(modelFlux[rowStart + j]);
        }

        var kernelWithStream = gpu.accelerator.LoadAutoGroupedKernel<Index1D, ArrayView<float>, ArrayView<float>, ArrayView<float>, ArrayView<float>, int>(CalcChiKernel);

        Vector output = new(gpu, new float[new_mod_flux.Count]);
        Vector new_mod_flux_vec = new(gpu, new_mod_flux.ToArray());

        using (GpuScope.Begin(output, new_mod_flux_vec, Flux, Error))
        {
            MemoryBuffer1D<float, Stride1D.Dense>
                buffer = output.GetBuffer(),
                buffer2 = new_mod_flux_vec.GetBuffer(),
                fluxBuf = Flux.GetBuffer(),
                errorBuf = Error.GetBuffer();

            kernelWithStream(gpu.accelerator.DefaultStream, buffer.IntExtent, buffer.View, buffer2.View, fluxBuf.View, errorBuf.View, model_len);
            gpu.accelerator.Synchronize();
        }

        ReadOnlySpan<float> chiFlat = output.RetrieveReadOnlySpan();
        float[] chis = new float[models];
        for (int i = 0; i < models; i++)
        {
            float sum = 0f;
            int rowStart = i * model_len;
            for (int j = 0; j < model_len; j++)
                sum += chiFlat[rowStart + j];
            chis[i] = sum;
        }

        return chis;
    }
    // Kernel
    static void CalcChiKernel(Index1D index, ArrayView<float> Output, ArrayView<float> InputMod_flux, ArrayView<float> InputDat_flux, ArrayView<float> InputDat_err, int model_len)
    {
        Output[index] = XMath.Pow((InputMod_flux[index] - InputDat_flux[index % model_len]) / InputDat_err[index % model_len], 2f);
    }
}
