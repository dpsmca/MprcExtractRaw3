namespace MprcExtractRaw3
{
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Collections.Specialized;
    using System.IO;
    using System.Linq;
    using System.Runtime.InteropServices;
    // using System.Runtime.Remoting.Messaging;
    using System.Text.RegularExpressions;
    using static System.Net.WebRequestMethods;

    // structure-like to store a peptide XIC
    public class PeptideChromatogram
    {
        // chrom = PeptideChromatogram(1.0, 4.0, 5.0, OrderedDict())
        // chrom.XICs[1] = [13., 14.]
        public double precursorMZ; // precursor m/z of the peptide
        // mz window for looking up the precursor m/z in the
        // spectral data. Window is defined as precursorMZ+/-ppmMassTol
        public double minMz;
        public double maxMz;
        /*
         * data storage of the XIC
         * data storage for (scan#, BasePeakXIC, and ticXIC). Storage is indexed by scan number
         * BasePeakXIC is the maximum intensity in the given range
         * ticXIC is the total intensity of all data points in the given range
         */
        //public Dictionary<int, Tuple<double, double>> XICs;  //  # PeakData type [scan_num] => (BasePeakXIC, ticXIC) i.e. first, second
        public OrderedDictionary XICs;  //  # PeakData type [scan_num] => (BasePeakXIC, ticXIC) i.e. first, second

        public PeptideChromatogram(double precursorMZ, double minMz, double maxMz, OrderedDictionary XICs)
        {
            this.precursorMZ = precursorMZ;
            this.minMz = minMz;
            this.maxMz = maxMz;
            this.XICs = XICs;
        }
    }

    public class ExtractorState
    {
        // Chromatograms
        public List<PeptideChromatogram> chromatograms;

        public ExtractorState(List<PeptideChromatogram> chromatograms)
        {
            this.chromatograms = chromatograms;
        }
    }

    // # Setup for extracting multiple chromatograms from the MS1 spectra
    public class ChromatogramExtractor
    {
        private string outputFileName;
        private string precursorMZsForXIC;
        private int massRounding;
        private float ppmMassTol;
        private ExtractorState state;
        private static Dictionary<string, int> headerDict = new Dictionary<string, int>()
        {
            { "Precursor m/z", 0 },
            { "m/z Window", 1 },
            { "Scan ID", 2 },
            { "BasePeakXIC", 3 },
            { "TICXIC", 4 },
            { "PeptideChromatogramNumber", 5 }  // only used for partial generated files
        };
        private static List<string> pdHeader = headerDict.Keys.Select(key => key.ToString()).ToList();

        public ChromatogramExtractor(string outputFileName, string precursorMZsForXIC, int massRounding, float ppmMassTol)
        {
            // Initialize the state. Called by the constructor.
            // File to save this output to
            this.outputFileName = outputFileName;
            this.precursorMZsForXIC = precursorMZsForXIC;  // a colon delimited list of precursor m/z values
            this.massRounding = massRounding;  // number of decimals for rounding the m/z values
            this.ppmMassTol = ppmMassTol; // mass tolerance for looking up the precursor m/z values in the spectral data
            // Internal state being built spectrum by spectrum
            this.state = new ExtractorState(new List<PeptideChromatogram>());

            // get the list of all precursor mzs. Each mz is separated by ':'
            // convert them into an array of doubles
            double[] precursors = new double[0];  // Initialize an empty array
            string[] precursorsList = precursorMZsForXIC.Split(':');
            foreach (string precursor in precursorsList)
            {
                if (precursor.Length >= 2)
                {
                    double precursorValue = double.Parse(precursor);
                    precursors = precursors.Append(precursorValue).ToArray();
                }
            }
            foreach (double precursorMZ in precursors)
            {
                double tol = (precursorMZ / 1.0e6) * this.ppmMassTol;
                double minMz = precursorMZ - tol;
                double maxMz = precursorMZ + tol;
                PeptideChromatogram chrom = new PeptideChromatogram(precursorMZ, minMz, maxMz, new OrderedDictionary());
                this.state.chromatograms.Add(chrom);
                // Console.WriteLine("INFO: #Precursor Peak Info " + Math.Round(precursorMZ, massRounding) + "\t" +
                //                  Math.Round(chrom.minMz, massRounding) + "-" + Math.Round(chrom.maxMz, massRounding));
            }
            // Console.WriteLine("INFO: Total number of precursors " + this.state.chromatograms.Count);
        }

        ~ChromatogramExtractor()
        {
            this.state.chromatograms.Clear();  // in case it's a big list
        }

        // UNIT TEST ME
        public string SavePartialRtcFile(bool toFile = true)
        {
            /*
                Save the extracted chromatogram data to a file or print to stdout.
                Output a text file containing the summed mass chromatogram
                :param to_file:  if True writes to self.outputFileName, else to stdout
                :return: string or an empty string (toFile == true)
                Designed to save only a small fraction of the total work (called in a separate thread).
             */
            List<string[]> pdRows = new List<string[]>();
            string pdStr = "";
            if (toFile)
            {
                // Console.WriteLine($"INFO: Writing peak list to {this.outputFileName}");
            }
            else
            {
                pdStr = string.Join("\t", pdHeader) + "\n";
            }

            int pcCntr = 0;
            foreach (PeptideChromatogram chromatogram in this.state.chromatograms)
            {
                int chromatogramLen = chromatogram.XICs.Count;
                if (chromatogramLen > 0)
                {
                    foreach (DictionaryEntry pdVal in chromatogram.XICs)
                    {
                        // C++ map<int, std::pair<double, double> >
                        int scanId = (int)pdVal.Key;
                        Tuple<double, double> value = (Tuple<double, double>)pdVal.Value;
                        double basePeakXIC = value.Item1;  // C++ pd->second.first
                        double ticXIC = value.Item2; // C++ pd->second.second
                        string[] pdRow = {
                            Math.Round(chromatogram.precursorMZ, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            Math.Round(chromatogram.minMz, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture) + "-" 
                            + Math.Round(chromatogram.maxMz, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            scanId.ToString(),
                            Math.Round(basePeakXIC, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            Math.Round(ticXIC, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            pcCntr.ToString()
                        };
                        if (toFile)
                        {
                            pdRows.Add(pdRow);
                        }
                        else
                        {
                            pdStr += string.Join("\t", pdRow) + "\n";
                        }
                    }
                }
                // We handle the chromatogramLen == 0 case when we concatenate the RTC files i.e. when we have the complete picture
                pcCntr++;
            }  //   foreachPeptideChromatogram chromatogram

            if (toFile)
            {
                using (StreamWriter writer = new StreamWriter(this.outputFileName))
                {
                    writer.WriteLine(string.Join("\t", pdHeader));
                    foreach (string[] row in pdRows)
                    {
                        writer.WriteLine(string.Join("\t", row));
                    }
                }
                return string.Empty;
            }
            else
            {
                return pdStr;
            }
        }

        // TODO: If the need for speed, then replace the following LINQ version with [].
        public void AddSpectrumOpt(int scanNumber, double[] mzs, double[] intensities)
        {
            /*
                Record a new spectrum
                :param scanNumber:
                :param mzs:
                :param intensities:
                :return:  fills in chromatogram.XICs[scanNumber][{0,1}]
                For optimization purposes, we assume that scan numbers are increasing.
             */
            // iterate over all candidate peptides
            foreach (PeptideChromatogram chromatogram in this.state.chromatograms)
            {
                // for each peak, bin it according to user input and keep a total of each bin
                var intensityCond = intensities.Select(x => x > 0).ToArray();
                
                // mz_cond = np.logical_and(chromatogram.minMz <= mzs, mzs <= chromatogram.maxMz)
                var mzCond = mzs.Select((mz, index) => chromatogram.minMz <= mz && mz <= chromatogram.maxMz).ToArray();
                
                var compoundCond = intensityCond.Zip(mzCond, (x, y) => x && y).ToArray();
                // handle boundary condition where we don't add to chromatogram.XICs[scanNumber]
                if (compoundCond.Any(b => b))
                {
                    // First == base peak, # Second == total ion current
                    // chromatogram.XICs[scanNumber] = [np.amax(intensities[compound_cond]), np.sum(intensities[compound_cond])]
                    Tuple<double, double> xics = new Tuple<double, double>(
                        intensities.Where((intensity, index) => compoundCond[index]).Max(),
                        intensities.Where((intensity, index) => compoundCond[index]).Sum()
                    );
                    chromatogram.XICs.Add(scanNumber, xics);
                }
            }
        }  // AddSpectrumOpt()

        public List<PeptideChromatogram> GetChromatograms()
        {
            return this.state.chromatograms;
        }

        public override string ToString()
        {
            return $"Chromatogram Extractor outputFileName = {this.outputFileName}, " +
                $"precursor m/z values = {this.precursorMZsForXIC}, mass rounding = {this.massRounding.ToString()}, ppm mass tolerance = {this.ppmMassTol.ToString()}{Environment.NewLine}" +
                $"{this.SavePartialRtcFile(false)}"; ;
        }

        public string ToDetailedString()
        {
            return $"Chromatogram Extractor({this.outputFileName}, " +
                $"{this.precursorMZsForXIC}, {this.massRounding.ToString()}, {this.ppmMassTol.ToString()}){Environment.NewLine}" + $"{this.SavePartialRtcFile(false)}";
        }

        public static bool SaveCompleteRtcFileWrapper(string output_dir, string rtc_base_filename, string rtc_output_filename,
           string FILE_EXT, string precursorMzs, int massRounding, float ppmMassTol, bool debug)
        {
            string[] rtc_files = Directory.GetFiles(output_dir, rtc_base_filename + ".*" + FILE_EXT);
            // Sort the filenames based on the extracted number part ex. Explorise240_with_uv.raw.rtc.1.tsv
            string pattern = @"(?<=raw\.rtc\.)\d+(?=\.tsv)";
            var sorted_rtc_files = rtc_files.OrderBy(filename =>
            {
                Match match = Regex.Match(filename, pattern);
                if (match.Success && int.TryParse(match.Value, out int number))
                {
                    return number;
                }
                return int.MaxValue; // Invalid numbers are treated as maximum value
            });
            if (sorted_rtc_files.Count() == 0)
            {
                Console.WriteLine("ERROR: No RTC files were generated for {0}", rtc_output_filename);
                return false;
            }

            ChromatogramExtractor chromatogramExtractor = new ChromatogramExtractor(rtc_output_filename, precursorMzs, massRounding, ppmMassTol);
            List<PeptideChromatogram> chromatograms = chromatogramExtractor.GetChromatograms();
            int num_chromatograms = chromatograms.Count();
            if (num_chromatograms == 0)
            {
                Console.WriteLine("ERROR: No chromatograms were configured for {0}", rtc_output_filename);
                return false;
            }    
            foreach (string file in sorted_rtc_files)
            {
                string[] lines = System.IO.File.ReadAllLines(file);
                lines = lines.Skip(1).ToArray();  // skip header
                for (int i = 0; i < lines.Length; i++)
                {
                    string line = lines[i];
                    string[] values = line.Split('\t');
                    if (values.Length != headerDict.Count)
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check columns", file);
                        return false;
                    }
                    int chromatogram_number = -1;
                    if (!Int32.TryParse(values[headerDict["PeptideChromatogramNumber"]], out chromatogram_number))
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check column 'PeptideChromatogramNumber'", file);
                        return false;
                    }
                    int scanNumber = -1;
                    if (!Int32.TryParse(values[headerDict["Scan ID"]], out scanNumber))
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check column 'Scan ID'", file);
                        return false;
                    }
                    double BasePeakXIC = 0.0;
                    if (!Double.TryParse(values[headerDict["BasePeakXIC"]], out BasePeakXIC))
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check column 'BasePeakXIC'", file);
                        return false;
                    }
                    double TICXIC = 0.0;
                    if (!Double.TryParse(values[headerDict["TICXIC"]], out TICXIC))
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check column 'TICXIC'", file);
                        return false;
                    }
                    if (chromatogram_number >= 0 && chromatogram_number < num_chromatograms)
                    {
                        Tuple<double, double> xics = new Tuple<double, double>(BasePeakXIC, TICXIC);
                        chromatograms[chromatogram_number].XICs.Add(scanNumber, xics);
                    }
                    else
                    {
                        Console.WriteLine("ERROR: Bad file format for {0}.  Check column 'PeptideChromatogramNumber'", file);
                        return false;
                    }
                }
            }  // foreach file

            // Write out complete RTC file in correct order            
            chromatogramExtractor.SaveCompleteRtcFile();

            // cleanup
            if (!debug)
            {
                foreach (string file in sorted_rtc_files)
                {
                    try
                    {
                        System.IO.File.Delete(file);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("ERROR: deleting the temporary file: {0}", file);
                        Console.WriteLine(ex.Message);
                    }
                }
            }  // if not debug
            return true;
        }  // SaveCompleteRtcFileWrapper()
        private void SaveCompleteRtcFile()
        {
            /*
                Save the extracted chromatogram data to a file.
                Output a text file containing the summed mass chromatogram
                :return: None
             */
            List<string[]> pdRows = new List<string[]>();
            Console.WriteLine($"INFO: Writing peak list to {this.outputFileName}");
            foreach (PeptideChromatogram chromatogram in this.state.chromatograms)
            {
                int chromatogramLen = chromatogram.XICs.Count;
                if (chromatogramLen > 0)
                {
                    foreach (DictionaryEntry pdVal in chromatogram.XICs)
                    {
                        // C++ map<int, std::pair<double, double> >
                        int scanId = (int)pdVal.Key;
                        Tuple<double, double> value = (Tuple<double, double>)pdVal.Value;
                        double basePeakXIC = value.Item1;  // C++ pd->second.first
                        double ticXIC = value.Item2; // C++ pd->second.second
                        string[] pdRow = {
                            Math.Round(chromatogram.precursorMZ, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            Math.Round(chromatogram.minMz, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture) +
                            "-" + Math.Round(chromatogram.maxMz, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            scanId.ToString(),
                            Math.Round(basePeakXIC, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                            Math.Round(ticXIC, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                        };
                        pdRows.Add(pdRow);
                    }
                }
                else 
                {
                    string[] pdRow = {
                        Math.Round(chromatogram.precursorMZ, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                        Math.Round(chromatogram.minMz, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture) + 
                            "-" + Math.Round(chromatogram.maxMz, massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                        "1",
                        Math.Round(0.00, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                        Math.Round(0.00, this.massRounding).ToString("0.00", System.Globalization.CultureInfo.InvariantCulture),
                    };
                    pdRows.Add(pdRow);
                }
            }  //   foreachPeptideChromatogram chromatogram
            using (StreamWriter writer = new StreamWriter(this.outputFileName))
            {
                List<string> pdHeaderFinal = new List<string>(pdHeader);
                pdHeaderFinal.Remove(pdHeaderFinal.Last());
                writer.WriteLine(string.Join("\t", pdHeaderFinal));
                foreach (string[] row in pdRows)
                {
                    writer.WriteLine(string.Join("\t", row));
                }
            }
        }  // SaveCompleteRtcFile()

    }  // ChromatogramExtractor

}  // namespace
