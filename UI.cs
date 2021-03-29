using System;
using System.Collections.Generic;
using System.Linq;


namespace MachineLearningSpectralFittingCode
{
    class UI
    {
        public UI(Config config)
        {
            RunUI(config);
        }

        private void RunUI(Config config)
        {
            bool running = true;
            Console.WriteLine("\tUserInterface 1.0 : ");
            while (running)
            {
                Console.WriteLine("\nPending Input ...");
                string input = Console.ReadLine();
                running = Commands(input, config);
            }
        }

        private bool Commands(string cmd, Config config)
        {
            switch (cmd.ToLower())
            {
                case "exit":
                    return false;

                case "help":
                    Help();
                    return true;
                case string s when cmd.Contains("-m"):
                    SetModel(s, config);
                    return true;



                default:
                    Console.WriteLine("\nWarning UNKNOWN input, either enter a valid input or " +
                        "type exit to close the program or type help for a detailed list of valid inputs.");
                    return true;
            }
        }


        private void Help()
        {
            Console.WriteLine("\n exit - will exit the user interface");
            Console.WriteLine("\n help - displays a list of valid user inputs");

            Console.WriteLine("\n -m/r - displays the current Model being used");
            Console.WriteLine("\n -m [model name] - allows modification of the current model used");


        }

        private void SetModel(string s, Config config)
        {
            byte[] modelkeys = new byte[6]
            {
                0b00001_001,
                0b00010_001,
                0b00100_001,
                0b01000_001,

                0b00001_010,
                0b00010_010,
            };
            string[] modelnames = new string[6]
            {
                "m11-miles",
                "m11-stelib",
                "m11-elodie",
                "m11-marcs",
                "mastar-th",
                "mastar-e",
            };

            if (s.Contains("-m/r"))
            {
                Console.WriteLine(
                    $"Current Model in use is : {modelnames[Array.IndexOf(modelkeys, config.Model_Key)]} with key : {config.Model_Key}");
                return;
            }

            string[] cmd = s.Split(" ");
            int idx = Array.IndexOf(cmd, "-m");

            // Call Config Class to Set model
            string model = "";
            try
            {
                model = cmd[idx + 1];
                config.Model_Key = modelkeys[Array.IndexOf(modelnames,model)];
                Console.WriteLine($"Current Model in use is {model} with key {config.Model_Key}");
            }
            catch (IndexOutOfRangeException)
            {
                Console.WriteLine("No Valid Input Provided For \"-m\" - Index Out Of Bounds ");
            }

            

        }



    }



}
