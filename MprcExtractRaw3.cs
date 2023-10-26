namespace MprcExtractRaw3
{
    using System;
    using System.Data;
    using System.Collections.Generic;
    using System.Collections.Specialized;
    using System.Diagnostics;
    using System.IO;
    using System.Linq;
    using System.Text.RegularExpressions;
    using System.Threading.Tasks;

    // These are the libraries necessary to read the Thermo RAW files.  The Interfaces
    // library contains the extension for accessing the scan averaging and background
    // subtraction methods.
    using ThermoFisher.CommonCore.Data;
    using ThermoFisher.CommonCore.Data.Business;
    using ThermoFisher.CommonCore.Data.FilterEnums;
    using ThermoFisher.CommonCore.Data.Interfaces;
    using ThermoFisher.CommonCore.RawFileReader;
    using System.Collections;
    using static global::MprcExtractRaw3.MprcUtils;
    using System.Security.Policy;
    using System.Collections.Concurrent;

    /// <summary>
    /// </summary>
    internal static class MprcExtractRaw3
    {
        const string FILE_EXT = ".tsv";

        private static float percentDone = 0.0F;
        private static object percentDoneLockObject = new object();

        /// <summary>
        /// </summary>
        /// <param name="args">
        /// The command line arguments for this program.  
        /// RAW filename            1st argument    (full path)
        /// Spectra filename        2nd             (full path, dir must exist)
        /// RTC filename            3rd             (full path, dir must exist)
        /// Chroma Gif filename     4th             (full path, dir must exist)
        /// Precursor Mzs           5th             string colon separated
        /// Mass Rounding            6th             int
        /// PPM Mass Tolerance      7th             float
        /// Debug flag              8th             (0 or 1=True)
        /// </param>
        /// WARNING!!! My RAW file handles should be closed prior to calling me from Python.
        private static void Main(string[] args)
        {
            // Get the memory used at the beginning of processing
            Process processBefore = Process.GetCurrentProcess();
            long memoryBefore = processBefore.PrivateMemorySize64 / 1024;
            var watch = System.Diagnostics.Stopwatch.StartNew();
            string paramsFile = string.Empty;
            try
            {
                string filename = string.Empty;
                string spectra_output_filename = string.Empty;
                string rtc_output_filename = string.Empty;
                string chroma_output_filename = string.Empty;
                string instrument_method_fileName = string.Empty;
                bool debug = false;
                // The following is an example and will be overwritten.
                string precursorMzs = "493.7683:613.3167:496.2867:451.2834:422.7363:695.8323:586.8003:773.8955:558.3259:801.4115:745.3924:498.8018:573.3025:680.3735:787.4212";
                int massRounding = 2;
                float ppmMassTol = 10.0F;
                int numRequestedCpus = 0;
                bool hasAllSpectraParms = true;
                if (args.Length == 2)
                {
                    string paramsFileKey = args[0];
                    if (!string.Equals(paramsFileKey, Config.PARAMS_FILE, StringComparison.InvariantCulture))
                    {
                        Console.WriteLine("ERROR: Only parameter {0} is supported right now", Config.PARAMS_FILE);
                        Environment.Exit(1);
                    }
                    paramsFile = args[1];
                    if (!File.Exists(paramsFile))
                    {
                        Console.WriteLine("ERROR: Parameters file does not exist, check path: " + paramsFile);
                        Environment.Exit(1);
                    }
                    OrderedDictionary parmsDict = GetParamsFromFile(paramsFile);
                    if (!HasOption(parmsDict, Config.DATA))
                    {
                        Console.WriteLine("ERROR: Parameters file does contain {0}", Config.DATA);
                        Environment.Exit(1);
                    }

                    filename = (string)GetOption(parmsDict, Config.RAW_FILE);
                    if (string.IsNullOrEmpty(filename))
                    {
                        Console.WriteLine("ERROR: No RAW file specified!");
                        Environment.Exit(1);
                    }
                    // Check to see if the specified RAW file exists
                    if (!File.Exists(filename))
                    {
                        Console.WriteLine(@"ERROR: The file doesn't exist in the specified location - " + filename);
                        Environment.Exit(1);
                    }
                    if (HasOption(parmsDict, Config.SPECTRA_FILE))
                    {
                        spectra_output_filename = (string)GetOption(parmsDict, Config.SPECTRA_FILE);
                    }
                    else
                    {
                        hasAllSpectraParms = false;
                    }
                    if (HasOption(parmsDict, Config.CHROMATOGRAM_FILE))
                    {
                        chroma_output_filename = (string)GetOption(parmsDict, Config.CHROMATOGRAM_FILE);
                    }
                    else
                    {
                        hasAllSpectraParms = false;
                    }
                    object numProcs = GetOption(parmsDict, Config.FIND_POLYMERS_NUM_PROCESSORS);
                    numRequestedCpus = (numProcs != null) ? Convert.ToInt32(numProcs) : 0;

                    object dbg = GetOption(parmsDict, Config.DEBUG);
                    if (dbg != null)
                    {
                        if(dbg.ToString() == "1")
                        {
                            debug = true;
                        }
                    }
                    if (HasOption(parmsDict, Config.RTC_FILE))
                    {
                        rtc_output_filename = (string)GetOption(parmsDict, Config.RTC_FILE);
                    }
                    else
                    {
                        hasAllSpectraParms = false;
                    }
                    if (!string.IsNullOrEmpty(rtc_output_filename))
                    {
                        precursorMzs = (string)GetOption(parmsDict, Config.RTC_PRECURSOR_MZS);
                        if (!string.IsNullOrEmpty(precursorMzs))
                        {
                            object mr = GetOption(parmsDict, Config.RTC_MASS_ROUNDING);
                            massRounding = (mr != null) ? Convert.ToInt32(mr) : -1;
                            if (massRounding < 1 || massRounding > 16) 
                            {
                                Console.WriteLine("ERROR: parameter {0} must be in the range [1, 16]", Config.RTC_MASS_ROUNDING);
                                Environment.Exit(1);
                            }
                            object ppmmt = GetOption(parmsDict, Config.RTC_PPM_MASS_TOL);
                            ppmMassTol = (ppmmt != null) ? Convert.ToSingle(ppmmt) : -1.0F;
                            if (ppmMassTol < 0)
                            {
                                Console.WriteLine("ERROR: parameter {0} must be must be >= 0", Config.RTC_PPM_MASS_TOL);
                                Environment.Exit(1);
                            }
                        }
                        else { hasAllSpectraParms = false; }
                    } // if rtc_output_filename
                    else { hasAllSpectraParms = false; }
                    instrument_method_fileName = (string)GetOption(parmsDict, Config.INSTRUMENT_METHOD_FILE);
                } // if args length
                else
                {
                    displayHelp();
                    Environment.Exit(1);
                }

                // ORIG single threaded example
                // Create the IRawDataPlus object for accessing the RAW file
                // var rawFile = RawFileReaderAdapter.FileFactory(filename);

                // MT interface
                IRawFileThreadManager myThreadManager = RawFileReaderFactory.CreateThreadManager(filename);
                IRawDataPlus rawFile = myThreadManager.CreateThreadAccessor();
                if (!rawFile.IsOpen || rawFile.IsError)
                {
                    Console.WriteLine("ERROR: Unable to access the RAW file using the RawFileReader class!");
                    Environment.Exit(1);
                }
                // Check for any errors in the RAW file
                if (rawFile.IsError)
                {
                    Console.WriteLine("ERROR: Error opening ({0}) - {1}", rawFile.FileError, filename);
                    Environment.Exit(1);
                }
                // Check if the RAW file is being acquired
                if (rawFile.InAcquisition)
                {
                    Console.WriteLine("ERROR: RAW file still being acquired - " + filename);
                    Environment.Exit(1);
                }

                // Get the number of instruments (controllers) present in the RAW file and set the 
                // selected instrument to the MS instrument, first instance of it
                Console.WriteLine("INFO:  The RAW file has data from {0} instruments", rawFile.InstrumentCount);

                // Used for general and instrument method calls
                rawFile.SelectInstrument(Device.MS, 1);
                if (!string.IsNullOrEmpty(instrument_method_fileName))
                {
                    rawFile.IncludeReferenceAndExceptionData = true;  // keep this visible for other calls to see
                    // NOTE-- UNCOMMENT when we can get Instrument Methods Data working under this project on Linux
                    // extractInstrumentMethodData(rawFile, instrument_method_fileName);
                }

                // Get the first and last scan from the RAW file
                int firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum;
                int lastScanNumber = rawFile.RunHeaderEx.LastSpectrum;
                int totalNumScans = lastScanNumber - firstScanNumber + 1;
                Console.WriteLine("INFO: scan numbers: {0}, {1}", firstScanNumber, lastScanNumber);
                rawFile.Dispose();   // Close (dispose) this thread's RAW file handle
                bool success = true;
                if (hasAllSpectraParms)
                {
                    success = extractPerSpectrumData(myThreadManager,
                        firstScanNumber, lastScanNumber, totalNumScans,
                        spectra_output_filename, rtc_output_filename, chroma_output_filename,
                        precursorMzs, massRounding, ppmMassTol,
                        numRequestedCpus, debug);
                }
                else
                {
                    Console.WriteLine("WARNING: Missing spectral parameters from {0}. Not extracting spectral data.", paramsFile);
                }
                // not sure if order matters here-- dispose the raw file first?
                myThreadManager.Dispose();
                Console.WriteLine();
                Console.WriteLine("INFO: Closing " + filename);
    
                watch.Stop();
                Console.WriteLine($"INFO: Execution Time: {watch.ElapsedMilliseconds/1000.0} sec");
                if (!success )
                {
                    Environment.Exit(1);
                }
            }  // try
            catch (FormatException ex)
            {
                Console.WriteLine("ERROR: {0}  Check {1}", ex.Message, paramsFile);
                Environment.Exit(1);
            }
            catch (Exception ex)
            {
                Console.WriteLine("ERROR: General exception - " + ex.Message);
                Environment.Exit(1);
            }

            // Get the memory used at the end of processing
            Process processAfter = Process.GetCurrentProcess();
            long memoryAfter = processAfter.PrivateMemorySize64 / 1024;

            Console.WriteLine();
            Console.WriteLine("INFO: Memory Usage:");
            Console.WriteLine("   Before {0} kb, After {1} kb, Extra {2} kb", memoryBefore, memoryAfter, memoryAfter - memoryBefore);
            Environment.Exit(0);
        }

        /// <summary>
      
        private static IEnumerable<int> SteppedIntegerList(int startIndex,  int endEndex, int stepSize)
        {
            for (int i = startIndex; i < endEndex; i += stepSize)
            {
                yield return i;
            }
        }  // SteppedIntegerList

        private static void displayHelp()
        {
            Console.WriteLine("Invalid arguments. Usage {0} (version {1}) --params <params_file>", 
                System.AppDomain.CurrentDomain.FriendlyName, Config.MPRC_EXTRACTRAW3_VERSION);
        }

        static bool extractPerSpectrumData(IRawFileThreadManager myThreadManager, 
            int first_scan_number, int last_scan_number, int total_num_scans,
            string spectra_output_filename, string rtc_output_filename, string chroma_output_filename,
            string precursorMzs, int massRounding, float ppmMassTol, 
            int numRequestedCpus, bool debug)
        {
            // get basename of raw file to derive 
            // 718213_Chup_20230421_Ex1_NDNF_S1.raw.spectra
            // 718213_Chup_20230421_Ex1_NDNF_S1.raw.rtc
            // 718213_Chup_20230421_Ex1_NDNF_S1.raw.chroma
            // NOTE-- Assume that we are writing to the same output dirtory for all related spectral files
            string output_dir = Directory.GetParent(spectra_output_filename).FullName + Path.DirectorySeparatorChar;
            // BUG FIX 9/12/2023
            if (!Directory.Exists(output_dir))
            {
                Console.WriteLine("ERROR: Unable to access directory {0}", output_dir);
                return false;
            }
            string spectra_base_filename = Path.GetFileNameWithoutExtension(spectra_output_filename);
            string rtc_base_filename = Path.GetFileNameWithoutExtension(rtc_output_filename);
            // Y:\data\Mayo\Proteo\new_huge_exploris\csharp\200ngHela_20230515_Ex1_1sAuto_SortHighestMz_Inj27.raw.chroma.{gif, png}
            string chroma_data_base_filename = Path.GetFileNameWithoutExtension(chroma_output_filename);

            // Extract low and high mass for first scan number (used for Chromatogram maps)
            IRawDataPlus rawFileUtil = myThreadManager.CreateThreadAccessor();
            IFileError fileErrorUtil = rawFileUtil.FileError;
            if (fileErrorUtil.HasError || fileErrorUtil.HasWarning)
            {
                Console.WriteLine("ERROR: Unable to access the Utility RAW file using the IRawDataPlus class! {0}, {1}",
                    fileErrorUtil.ErrorMessage, fileErrorUtil.WarningMessage);
                return false;
            }
            rawFileUtil.SelectInstrument(Device.MS, 1);
            (double lowMass, double highMass) = MprcUtils.GetLowHighMassFromFirstScanNumber(rawFileUtil, first_scan_number);
            rawFileUtil.Dispose();

            int num_processors = numRequestedCpus;
            int num_logical_processors = Environment.ProcessorCount;
            if (num_processors == 0)
            {
                num_processors = num_logical_processors;
                Console.WriteLine("INFO: Number Of Logical Processors: {0}", num_processors);
                num_processors /= 4;  // Using a quarter of the available cores produced the fastest result for two 122K scans.
            }
            num_processors = Math.Min(num_processors, num_logical_processors - 1);  // save room for possible extra work
            // ex.)  11 scans using 2 processors => groupSize=5 remainSize=1 [1,10] [11]
            // ex.)  15 scans using 4 processors => groupSize=3 remainSize=3 [1,12] [13-15]  
            int groupSize = total_num_scans / num_processors;
            int remainSize = total_num_scans % num_processors;
            // This will work bec. we use less than the totoal number of processors
            if (remainSize > 0)
            {
                num_processors++;
            }
            Console.WriteLine("INFO: Using Max Degree Of Parallelism {0}", num_processors);
            var options = new ParallelOptions()
            {
                MaxDegreeOfParallelism = num_processors
            };
            // ConcurrentBag to store the results safely
            ConcurrentBag<double> maxLogIntensities = new ConcurrentBag<double>();
            Parallel.ForEach(SteppedIntegerList(first_scan_number, last_scan_number + 1, groupSize), options, scan_num =>
            {
                IRawDataPlus rawFileMT = myThreadManager.CreateThreadAccessor();
                IFileError fileErrorMT = rawFileMT.FileError;
                if (fileErrorMT.HasError || fileErrorMT.HasWarning)
                {
                    Console.WriteLine("ERROR: Unable to access the RAW file using the IRawDataPlus class! {0}, {1}",
                        fileErrorMT.ErrorMessage, fileErrorMT.WarningMessage);
                    return;  // cannot return a boolean here, hopefully ok
                }
                bool firstSpectrum = false;
                if (scan_num == first_scan_number)
                    firstSpectrum = true;

                int actualGroupSize = groupSize;
                if (scan_num == total_num_scans - remainSize + 1)
                {
                    actualGroupSize = remainSize;
                }
                double maxLogIntensity = extractPerSpectrumDataTask(rawFileMT, scan_num, actualGroupSize, total_num_scans,
                    spectra_base_filename, rtc_base_filename, chroma_data_base_filename, output_dir, firstSpectrum,
                    precursorMzs, massRounding, ppmMassTol, lowMass, highMass);
                maxLogIntensities.Add(maxLogIntensity);
            });  // parfor

            // Thread join step to concatenate the separate thread files.
            ConcatenateSpectralThreadFiles(output_dir, spectra_base_filename, spectra_output_filename, debug);
            bool success1 = ChromatogramExtractor.SaveCompleteRtcFileWrapper(output_dir, rtc_base_filename, rtc_output_filename,
                FILE_EXT, precursorMzs, massRounding, ppmMassTol, debug);

            // Find the maximum value among the results
            double theMaxLogIntensity = maxLogIntensities.Max();
            bool success2 = ChromatogramMap.DumpEqualizedWrapper(output_dir, chroma_data_base_filename, chroma_output_filename, FILE_EXT,
                Config.CHROMATOGRAM_MZ_BINS, lowMass, highMass, theMaxLogIntensity, debug);
            return (success1 && success2);
        }

        // I calculate a range of scan numbers in parallel.
        static double extractPerSpectrumDataTask(IRawDataPlus raw_file, int first_scan_number, int step_scan_num, int total_num_scans,
            string spectra_base_filename, string rtc_base_filename, string chroma_data_base_filename, string output_dir, bool firstSpectrum,
              string precursorMzs, int massRounding, float ppmMassTol, double firstLowMass, double firstHighMass)
        {
            double[] mzs = new double[0];
            double[] intensities = new double[0];
            // The following can be created in each task (thread) for its unit of work.
            PolymerDetection pd = new PolymerDetection();
           
            int last_scan_number = first_scan_number + step_scan_num;  // upper bounded
            //Console.WriteLine("INFO: START In new thread processing scan #s: {0}, {1}, {2}", 
            //    first_scan_number, last_scan_number - 1, firstSpectrum.ToString());
            // 718213_Chup_20230421_Ex1_NDNF_S1.raw.spectra
            // DESIGN-- we could also store in a subdir, don't forget cleanup if not debug
            // ex. 718213_Chup_20230421_Ex1_NDNF_S1.raw.spectra.1.tsv
            string my_spectra_filename = output_dir + spectra_base_filename + "." + first_scan_number.ToString() + FILE_EXT;
            string my_rtc_filename = output_dir + rtc_base_filename + "." + first_scan_number.ToString() + FILE_EXT;
            string my_chroma_data_filename = output_dir + chroma_data_base_filename + "." + first_scan_number.ToString() + FILE_EXT;

            ChromatogramExtractor chromatogramExtractor = new ChromatogramExtractor(my_rtc_filename, precursorMzs, massRounding, ppmMassTol);

            List<string> spectrumHeader = new List<string> { Config.SCAN_ID, Config.PARENT_MZ, Config.TOTAL_ION_CURRENT, Config.RETENTION_TIME, Config.MS_LEVEL,
                       Config.PARENT_SCAN, Config.CHILD_SCANS, Config.ION_INJECTION_TIME_MS, Config.CYCLE_TIME_SECONDS, Config.ELAPSED_TIME_SECONDS,
                       Config.DEAD_TIME_SECONDS, Config.TIME_TO_NEXT_SCAN_SECONDS, Config.LOCK_MASS_FOUND, Config.LOCK_MASS_SHIFT, Config.CON_I,
                       Config.CON_A, Config.CON_B, Config.CON_C, Config.CON_D, Config.CON_E,
                       Config.DISSOCIATION_TYPE, Config.POLYMER_SEGMENT, Config.POLYMER_OFFSET, Config.POLYMER_SCORE, Config.POLYMER_PVALUE };
            foreach (int mass in Config.POLYMER_MASSES)
            {
                spectrumHeader.Add(string.Format(Config.POLYMER_FOR_MASS_OFFSET, mass.ToString()));
                spectrumHeader.Add(string.Format(Config.POLYMER_FOR_MASS_SCORE, mass.ToString()));
                spectrumHeader.Add(string.Format(Config.POLYMER_FOR_MASS_PVALUE, mass.ToString()));
            }
            List<string> spectrumHeaderAdd = new List<string>
            {
                Config.BASE_PEAK_MZ, Config.BASE_PEAK_INTENSITY, Config.SECOND_PEAK_MZ, Config.SECOND_PEAK_INTENSITY,
                Config.SOURCE_CURRENT, Config.VACUUM_ION_GAUGE, Config.VACUUM_CONVECTRON_GAUGE,
                Config.FT_VACUUM_PENNING_GAUGE, Config.FT_VACUUM_PIRANI_GAUGE1, Config.ION_MULTIPLIER1,
                Config.ION_MULTIPLIER2, Config.FT_CE_MEASURE_VOLTAGE, Config.FT_ANALYZER_TEMP
            };
            spectrumHeader.AddRange(spectrumHeaderAdd);
            ChromatogramMap map = null;  // ChromatogramMap
            List<List<string>> spectrumDataRows = new List<List<string>>();  // may need mixed types in this list
            raw_file.SelectInstrument(Device.MS, 1);
            int num_empty_scans = 0;
            
            for (int scan_num = first_scan_number; scan_num < last_scan_number; scan_num++)
            {
                // init local vars
                double parentMz = 0.0;
                double tic = 0.0;
                double retentionTime = 0.0;
                int msLevel = 1;
                int parentScan = 0;
                double lowMass = 0.0;  // currently unused as we precalc it
                double highMass = 0.0;  // currently unused as we precalc it
                int childScans = 0;
                double ionInjectionTimeMs = 0.0;
                double cycleTimeSeconds = 0.0;
                double elapsedTimeSeconds = 0.0;
                double deadTimeSeconds = 0.0;
                double timeToNextScanSeconds = 0.0;
                bool lockMassFound = false;
                double lockMassShift = 0.0;
                double conI = 0.0;
                double conA = 0.0;
                double conB = 0.0;
                double conC = 0.0;
                double conD = 0.0;
                double conE = 0.0;
                int polymerSegment = -1;
                int polymerOffset = -1;
                double polymerScore = 0.0;

                double polymerPValue = 1.0;
                double basePeakMz = 0.0;
                double basePeakIntensity = 0.0;
                double secondPeakMz = 0.0;
                double secondPeakIntensity = 0.0;
                double sourceCurrent = 0.0;
                double vacuumIonGauge = 0.0;
                double vacuumConvectronGauge = 0.0;
                double ftVacuumPenningGauge = 0.0;
                double ftVacuumPiraniGauge1 = 0.0;
                double ionMultiplier1 = 0.0;
                double ionMultiplier2 = 0.0;
                double ftCeMeasureVoltage = 0.0;
                double ftAnalyzerTemp = 0.0;

                (tic, retentionTime, lowMass, highMass, childScans, ionInjectionTimeMs,
                 cycleTimeSeconds, elapsedTimeSeconds,
                 timeToNextScanSeconds,
                 lockMassFound, lockMassShift,
                 conI, conA, conB, conC, conD, conE,
                 sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge, ftVacuumPiraniGauge1,
                 ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp, basePeakMz, basePeakIntensity) = MprcUtils.GetScanInfo(raw_file, scan_num);

                msLevel = MprcUtils.GetMSLevel(raw_file, scan_num);
                if (msLevel != 1)
                {
                    parentMz = MprcUtils.GetParentMz(raw_file, scan_num);
                }
                else
                {
                    parentMz = 0.0;
                }

                if (map == null && !string.IsNullOrEmpty(my_chroma_data_filename))
                {
                    map = new ChromatogramMap(Config.CHROMATOGRAM_MZ_BINS, firstLowMass, firstHighMass);
                }
                bool calculatePolymers = false;
                double mz = 0.0;
                int charge = 0;
                msLevel = MprcUtils.GetMSLevel(raw_file, scan_num);
                if (msLevel == 1)
                {
                    deadTimeSeconds = cycleTimeSeconds - elapsedTimeSeconds;
                    parentScan = 0;
                    bool success = false;
                    if (map != null || chromatogramExtractor != null)
                    {
                        (success, mzs, intensities) = MprcUtils.GetRawData(raw_file, scan_num);
                        if (success && map != null)
                        {
                            map.AddSpectrum(mzs, intensities);
                        }
                    }
                    
                    if (success && chromatogramExtractor != null) 
                    {
                        chromatogramExtractor.AddSpectrumOpt(scan_num, mzs, intensities);
                    }

                    if (!success)
                    {
                        num_empty_scans += 1;
                    }
                }
                else
                {
                    if (!string.IsNullOrEmpty(my_spectra_filename))
                    {
                        deadTimeSeconds = timeToNextScanSeconds - elapsedTimeSeconds;
                        parentScan = MprcUtils.GetParentScan(raw_file, scan_num);
                        bool success = false;
                        (success, mzs, intensities, mz, charge) = MprcUtils.GetRawDataWrapper(raw_file, scan_num);
                        calculatePolymers = success;
                        if(success)
                        {
                            (polymerSegment, polymerOffset, polymerScore, polymerPValue) = PolymerScore.GetPolymerScore(mzs, intensities, mz, charge, 
                                Config.MIN_SEGMENT_SIZE, Config.MAX_SEGMENT_SIZE, pd);
                            (secondPeakMz, secondPeakIntensity) = MprcUtils.GetBasePeak(mzs, intensities, basePeakMz, basePeakIntensity, Config.SECOND_PEAK_MIN_DISTANCE_FROM_BASE);
                        }
                        else
                        {
                            num_empty_scans += 1;
                        }
                    }
                }  // else not ms level 1

                if (!string.IsNullOrEmpty(my_spectra_filename))
                {
                    // not sure if we should cast all of these data items to strings wrt decimal places?
                    List<string> spectrumDataRow = new List<string>
                    {
                        scan_num.ToString(),
                        parentMz != 0.0 ? parentMz.ToString("F8") : "",
                        tic.ToString("F8"),
                        retentionTime.ToString("F8"),
                        msLevel.ToString(),
                        parentScan != 0 ? parentScan.ToString() : "",
                        childScans.ToString(),
                        firstSpectrum && ionInjectionTimeMs == (int)ionInjectionTimeMs ?
                            ionInjectionTimeMs.ToString("F1") : ionInjectionTimeMs.ToString("F8"),
                        cycleTimeSeconds.ToString("F8"),
                        elapsedTimeSeconds.ToString("F8"),
                        deadTimeSeconds.ToString("F8"),
                        timeToNextScanSeconds.ToString("F8"),
                        lockMassFound ? "1" : "0",
                        lockMassShift.ToString("F8"),
                        conI.ToString("F8"),
                        conA.ToString("F8"),
                        conB.ToString("F8"),
                        conC.ToString("F8"),
                        conD.ToString("F8"),
                        conE.ToString("F8"),
                        MprcUtils.GetDissociationType(raw_file, scan_num),
                        polymerSegment.ToString(),
                        polymerOffset.ToString(),
                        polymerScore.ToString("F8"),
                        firstSpectrum && polymerPValue == (int)polymerPValue ?
                            polymerPValue.ToString("F1") : polymerPValue.ToString("F8")
                    };

                    List<int> POLYMER_MASSES = new List<int> { 162 };
                    for (int polymerId = 0; polymerId < POLYMER_MASSES.Count; polymerId++)
                    {
                        int mass = POLYMER_MASSES[polymerId];
                        int polymerForMassOffset = -1;
                        double polymerForMassScore = 0;
                        double polymerForMassPValue = 1.0;
                        if (calculatePolymers)
                        {
                            int segment = POLYMER_MASSES[polymerId];
                            int polymerForMassSegment = -1;
                            (polymerForMassSegment, polymerForMassOffset, polymerForMassScore, polymerForMassPValue) = PolymerScore.GetPolymerScore(mzs, intensities, mz, charge,
                                segment, segment, pd);
                        }
                        spectrumDataRow.Add(polymerForMassOffset.ToString());
                        spectrumDataRow.Add(polymerForMassScore.ToString("F8"));
                        spectrumDataRow.Add( firstSpectrum && polymerForMassPValue == (int)polymerForMassPValue ?
                            polymerForMassPValue.ToString("F1") : polymerForMassPValue.ToString("F8") );
                    } // for polymerId
                    // Base peak and second most intense
                    if (firstSpectrum && Math.Abs(basePeakMz - 0.0) < double.Epsilon)
                    {
                        spectrumDataRow.Add(basePeakMz.ToString("0.0"));
                    }
                    else
                    {
                        spectrumDataRow.Add(basePeakMz.ToString("F8"));
                    }
                    if (firstSpectrum && Math.Abs(basePeakIntensity - 0.0) < double.Epsilon)
                    {
                        spectrumDataRow.Add(basePeakIntensity.ToString("0.0"));
                    }
                    else
                    {
                        spectrumDataRow.Add(basePeakIntensity.ToString("F8"));
                    }
                    if (firstSpectrum && Math.Abs(secondPeakMz - 0.0) < double.Epsilon)
                    {
                        spectrumDataRow.Add(secondPeakMz.ToString("0.0"));
                    }
                    else
                    {
                        spectrumDataRow.Add(secondPeakMz.ToString("F8"));
                    }
                    if (firstSpectrum && Math.Abs(secondPeakIntensity - 0.0) < double.Epsilon)
                    {
                        spectrumDataRow.Add(secondPeakIntensity.ToString("0.0"));
                    }
                    else
                    {
                        spectrumDataRow.Add(secondPeakIntensity.ToString("F8"));
                    }

                    // Status log
                    spectrumDataRow.Add(sourceCurrent.ToString("F8"));
                    spectrumDataRow.Add(vacuumIonGauge.ToString("F8"));
                    spectrumDataRow.Add(vacuumConvectronGauge.ToString("F8"));
                    spectrumDataRow.Add(ftVacuumPenningGauge.ToString("F8"));
                    spectrumDataRow.Add(ftVacuumPiraniGauge1.ToString("F8"));
                    spectrumDataRow.Add(ionMultiplier1.ToString("F8"));
                    spectrumDataRow.Add(ionMultiplier2.ToString("F8"));
                    spectrumDataRow.Add(ftCeMeasureVoltage.ToString("F8"));
                    spectrumDataRow.Add(ftAnalyzerTemp.ToString("F8"));

                    spectrumDataRows.Add(spectrumDataRow);
                }  // if my_spectra_filename

                if (firstSpectrum == true)
                    firstSpectrum = false;
            }  // for scan_num
           
            if (!string.IsNullOrEmpty(my_spectra_filename))
            {
                DataTable dataTable = new DataTable();
                // Create the columns in the DataTable based on the first row of the spectrumDataRows list
                for (int i = 0; i < spectrumHeader.Count; i++)
                {
                    dataTable.Columns.Add(spectrumHeader[i]); // Replace with your desired column names
                }
                // Add the data from spectrumDataRows to the DataTable
                foreach (List<string> rowData in spectrumDataRows)
                {
                    DataRow dataRow = dataTable.NewRow();
                    for (int i = 0; i < rowData.Count; i++)
                    {
                        dataRow[i] = rowData[i];
                    }
                    dataTable.Rows.Add(dataRow);
                }
                using (StreamWriter writer = new StreamWriter(my_spectra_filename))
                {
                    writer.WriteLine(string.Join("\t", dataTable.Columns.Cast<DataColumn>().Select(c => c.ColumnName)));
                    foreach (DataRow row in dataTable.Rows)
                    {
                        writer.WriteLine(string.Join("\t", row.ItemArray));
                    }
                }
            }  // if my_spectra_filename

            if (chromatogramExtractor != null)
            {
                chromatogramExtractor.SavePartialRtcFile();
            }
            if (map != null) 
            {
                map.SavePartialChromaMapFile(my_chroma_data_filename);
            }

            if (num_empty_scans != 0)
            { 
                Console.WriteLine("INFO: Number of empty scans = {0}", num_empty_scans);
            }
            raw_file.Dispose();

            float percentDoneSpectral = (float)step_scan_num / (float)total_num_scans * 100.0F;
            float progress = 0.0F;
            lock (percentDoneLockObject)
            {
                percentDone += percentDoneSpectral;
                progress = percentDone;
            }
            Console.WriteLine("INFO: processed scan #s: {0}, {1}, percent done={2}",
              first_scan_number, last_scan_number - 1, progress.ToString("0.00", System.Globalization.CultureInfo.InvariantCulture));

            if (map != null)
            {
                return map.GetMaxLogIntensity();
            }
            else { return 0; }
        }  // extractPerSpectrumDataTask()

        static void ConcatenateSpectralThreadFiles(string output_dir, string spectra_base_filename, string spectra_output_filename,
            bool debug)
        {
            string[] spectral_files = Directory.GetFiles(output_dir, spectra_base_filename + ".*" + FILE_EXT);
            // Sort the filenames based on the extracted number part ex. Explorise240_with_uv.raw.spectra.1.tsv
            string pattern = @"(?<=raw\.spectra\.)\d+(?=\.tsv)";
            var sorted_spectral_files = spectral_files.OrderBy(filename =>
            {
                Match match = Regex.Match(filename, pattern);
                if (match.Success && int.TryParse(match.Value, out int number))
                {
                    return number;
                }
                return int.MaxValue; // Invalid numbers are treated as maximum value
            });
            bool isFirstFile = true;
            using (StreamWriter writer = new StreamWriter(spectra_output_filename))
            {
                foreach (string file in sorted_spectral_files)
                {
                    // Read the contents of each file
                    string[] lines = File.ReadAllLines(file);

                    // Write the lines to the output file, excluding the header line if it's not the first file
                    if (!isFirstFile)
                    {
                        lines = lines.Skip(1).ToArray();
                    }
                    writer.WriteLine(string.Join(Environment.NewLine, lines));

                    // Set isFirstFile to false after processing the first file
                    if (isFirstFile)
                    {
                        isFirstFile = false;
                    }
                }  // foreach file
            } // using writer
            // cleanup
            if (!debug)
            {
                foreach (string file in sorted_spectral_files)
                {
                    try
                    {
                        File.Delete(file);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("ERROR: deleting the temporary file: {0}", file);
                        Console.WriteLine(ex.Message);
                    }
                }
            }  // if not debug
        }  // ConcatenateSpectralThreadFiles()

        public static OrderedDictionary GetParamsFromFile(string inputParamsFile)
        {
            /*
                :param inputParamsFile:
                :return: Ordered dictionary of parameter names => values
            */
            OrderedDictionary parmsDict = new OrderedDictionary();
            Console.WriteLine("Reading parameter file " + inputParamsFile + ".");
            string[] parmsList = File.ReadAllLines(inputParamsFile);
            foreach (string eachLine in parmsList)
            {
                string line = eachLine.Trim();  // strip white space
                if (line.StartsWith("#") || line == "")
                    continue;
                line = line.Trim('"');  // strip quotes off
                //parmsDict.Add(line, null);  // unsure why G added this?
            }
            int parmsListLen = parmsList.Length;
            for (int index = 0; index < parmsListLen; index++)
            {
                string item = parmsList[index];
                if (item == "--data")  // data has no subsequent value
                {
                    parmsDict[item] = null;
                }
                else if (index <= parmsListLen - 2 && item.StartsWith("--"))
                {
                    // If we are not at the end of the list, and we have a parameter (key)
                    // the next item is assumed to be the value (from a separate line in the file)
                    parmsDict[item] = parmsList[index + 1];
                }
            }
            foreach (DictionaryEntry entry in parmsDict)
            {
                string k = (string)entry.Key;
                object v = entry.Value;
                Console.WriteLine($"{k} = {v}");
            }
            return parmsDict;
        }

        public static bool HasOption(OrderedDictionary parmsDict, string option)
        {
            return parmsDict.Contains(option);
        }

        public static object GetOption(OrderedDictionary parmsDict, string option)
        {
            return parmsDict[option];  // <-- No default specified -- defaults to null
        }

        //public static void Main()
        //{
        //    string inputParamsFile = "path/to/params/file.txt";
        //    OrderedDictionary parmsDict = GetParamsFromFile(inputParamsFile);

        //    string option = "--data";
        //    bool hasOption = HasOption(parmsDict, option);
         

        //    object optionValue = GetOption(parmsDict, option);
           
        //}
    
    }  // class MprcExtractRaw3
}  // namespace
