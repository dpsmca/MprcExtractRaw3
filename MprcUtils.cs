namespace MprcExtractRaw3
{
    using System;
    using System.Collections.Generic;
    using System.Collections.Specialized;
    using System.Linq;
    using System.IO;

    using ThermoFisher.CommonCore.Data;
    using ThermoFisher.CommonCore.Data.Business;
    using ThermoFisher.CommonCore.Data.FilterEnums;
    using ThermoFisher.CommonCore.Data.Interfaces;
    using ThermoFisher.CommonCore.MassPrecisionEstimator;
    using ThermoFisher.CommonCore.RawFileReader;
    using System.Security.Cryptography;

    public static class MprcUtils
    {

        public class ScanTrailer
        {
            private Dictionary<string, string> _data;

            public int Length => _data.Count;

            public List<string> Labels => new List<string>(_data.Keys);

            public List<string> Values => new List<string>(_data.Values);

            public bool? AsBool(string key)
            {
                if (_data.ContainsKey(key))
                {
                    string strValue = _data[key].ToLower();

                    if (strValue == "on" || strValue == "true" || strValue == "yes")
                    {
                        return true;
                    }
                    return false;
                }
                return null;
            }

            public double? AsDouble(string key)
            {
                if (_data.ContainsKey(key))
                {
                    if (double.TryParse(_data[key], out double result))
                    {
                        return result;
                    }
                }
                return null;
            }

            public int? AsInt(string key)
            {
                if (_data.ContainsKey(key))
                {
                    if (int.TryParse(_data[key], out int result))
                    {
                        return result;
                    }
                }
                return null;
            }

            public int? AsPositiveInt(string key)
            {
                if (_data.ContainsKey(key))
                {
                    if (int.TryParse(_data[key], out int result))
                    {
                        return result >= 0 ? result : (int?)null;
                    }
                }
                return null;
            }

            public string AsString(string key)
            {
                return Get(key);
            }

            public string Get(string key)
            {
                if (_data.ContainsKey(key))
                {
                    return _data[key];
                }
                return null;
            }

            public bool Has(string key)
            {
                return _data.ContainsKey(key);
            }

            public ScanTrailer(LogEntry trailerData)
            {
                _data = new Dictionary<string, string>();
                for (int i = 0; i < trailerData.Length; i++)
                {
                    _data.Add(trailerData.Labels[i], trailerData.Values[i]);
                }
            }
    }

        public class MetadataWriter
        {
            private string _metadataFileName;

            /// <summary>
            /// Constructor.
            /// </summary>
            /// <param name="metadataFileName"></param>
            public MetadataWriter(string metadataFileName)
            {
                _metadataFileName = metadataFileName;
            }

            /// <summary>
            /// Write the RAW file metadata to file.
            /// <param name="rawFile">the RAW file object</param>
            /// </summary>
            public void WriteMetadata(IRawDataPlus rawFile)
            {
                WriteTextMetadata(rawFile);
            }

            /// <summary>
            /// Write the RAW file metadata to file.
            /// <param name="rawFile">the RAW file object</param>
            /// <param name="firstScanNumber">the first scan number</param>
            /// <param name="lastScanNumber">the last scan number</param>
            /// </summary>
            private void WriteTextMetadata(IRawDataPlus rawFile)
            {
                // File Properties
                var output = new List<string>();
                // Instrument method outputs in slightly different order.
                // May need to change encoding for Quameter-- https://docs.python.org/3/library/codecs.html#standard-encodings
                for (int i = 0; i < rawFile.InstrumentMethodsCount; i++)
                {
                    var methodText = rawFile.GetInstrumentMethod(i);
                    if (!String.IsNullOrEmpty(methodText))
                    {
                        var splitMethod = methodText.Split(new[] { "\n" }, StringSplitOptions.None);
                        foreach (var line in splitMethod)
                        {
                            var lineCleaned = line.Replace("\r", "");
                            output.Add(lineCleaned);
                        }
                    }
                }
                // Write the meta data to file
                File.WriteAllLines(_metadataFileName, output.ToArray());
            }
        }  // class MetadataWriter

        public static (bool success, double[] mzs, double[] intensities) GetRawData(IRawDataPlus rawFile, int scanNum)
        {
            int nArraySize = 0;
            rawFile.SelectInstrument(Device.MS, 1);

            // "You can use any of these methods (GetCentroidStream, GetSegmentedScan, Scan.FromFile)."-- Jim Shofstahl from Thermo Fisher                          
            ScanStatistics scanStatistics = rawFile.GetScanStatsForScanNumber(scanNum);
            double highMass = scanStatistics.HighMass;
            double[] mzs = new double[0];
            double[] intensities = new double[0];
            // the following is used to prevent duplicate code                                                                                                      
            double[] mzsCommon = new double[0];
            double[] intensitiesCommon = new double[0];

            // Check to see if the scan has centroid data or profile data.  Depending upon the                                                                      
            // type of data, different methods will be used to read the data.  While the ReadAllSpectra                                                             
            // method demonstrates reading the data using the Scan.FromFile method, generating the                                                                  
            // Scan object takes more time and memory to do, so that method isn't optimum.                                                                          
            if (scanStatistics.IsCentroidScan && scanStatistics.SpectrumPacketType == SpectrumPacketType.FtCentroid)
            {
                // Get the centroid (label) data from the RAW file for this scan                                                                                    
                CentroidStream centroidStream = rawFile.GetCentroidStream(scanNum, false);
                nArraySize = centroidStream.Length; // size of mass list array                                                                                     
                if (nArraySize == 0)
                {
                    return (false, mzs, intensities);
                }
                // Enum the spectral data (mass, intensity, charge values).                                                                                         
                // We needed a list() to convert from .NET []                                                                                                       
                mzsCommon = centroidStream.Masses;
                intensitiesCommon = centroidStream.Intensities;
                // we also have centroidStream.Charges[i]                                                                                                          
            }
            else
            {
                // Get the segmented (low res and profile) scan data                                                                                                
                SegmentedScan segmentScan = rawFile.GetSegmentedScanFromScanNumber(scanNum, scanStatistics);
                nArraySize = segmentScan.Positions.Length; // size of mass list array                                                                              
                if (nArraySize == 0)
                {
                    return (false, mzs, intensities);
                }
                mzsCommon = segmentScan.Positions;
                intensitiesCommon = segmentScan.Intensities;
            }

            // common logic                                                                                                                                         
            int[] idxsLeHighMass = Array.FindAll(
                Enumerable.Range(0, mzsCommon.Length).ToArray(),
                i => mzsCommon[i] <= highMass);
            mzs = Array.ConvertAll(idxsLeHighMass, i => mzsCommon[i]);
            intensities = Array.ConvertAll(idxsLeHighMass, i => intensitiesCommon[i]);
            if (mzs.Length == 0)
            {
                return (false, mzs, intensities);
            }
            else
            {
                return (true, mzs, intensities);
            }
        } // GetRawData

        public static (bool success, double[] mzs, double[] intensities, double mz, int charge) 
            GetRawDataWrapper(IRawDataPlus raw_file, int scan_num)
        {
            /*
            :param raw_file:
            :param scan_num:
            :return: success, mzs, intensities, mz, charge
            */
            double[] mzs = new double[0];
            double[] intensities = new double[0];
            bool success;
            (success, mzs, intensities) = GetRawData(raw_file, scan_num);  // should check return code
            (double mz, int charge) = GetMonoMZChargeFromHeader(raw_file, scan_num);
            return (success, mzs, intensities, mz, charge);
        }

        public static (double, int) GetMonoMZChargeFromHeader(IRawDataPlus rawFile, int scanNum)
        {
            double monoMZVal = 0.0;
            int csVal = 0;
            // May return null
            ScanTrailer trailerData = new ScanTrailer(rawFile.GetTrailerExtraInformation(scanNum));
            double ?monoMZ = trailerData.AsDouble("Monoisotopic M/Z:");
            // duplicated behavior of C++ m_xraw2_class->GetTrailerExtraValueForScanNum()
            if (monoMZ != null)
            {
                monoMZVal = (double)monoMZ;
            }
            int? cs = trailerData.AsPositiveInt("Charge State:");
            if (cs != null)
            {
                csVal = (int)cs;
            }
            return (monoMZVal, csVal);
        }

        public static double GetParentMz(IRawDataPlus raw_file, int scan_num)
        {
            double parent_mz = 0.0;
            string ch_mz = "";
            ScanStatistics scan_statistics = raw_file.GetScanStatsForScanNumber(scan_num);
            string ch_filter = scan_statistics.ScanType;
            int ms_level = GetMSLevel(raw_file, scan_num);
            int parent_count = 0;

            if (ms_level == 2)
            {
                int chNum = ch_filter.IndexOf('2');
                if (chNum != -1)
                {
                    chNum += 1;
                    while (chNum < ch_filter.Length && ch_filter[chNum] != '@')
                    {
                        ch_mz += ch_filter[chNum];
                        chNum += 1;
                    }
                }
            }
            else
            {
                int chNum = 0;
                for (int i = 0; i < ch_filter.Length; i++)
                {
                    if (ch_filter[i] == '@')
                    {
                        parent_count += 1;
                        if (parent_count <= (ms_level - 1))
                        {
                            chNum = i;
                            break;
                        }
                    }
                }
                while (chNum < ch_filter.Length && ch_filter[chNum] != ' ')
                {
                    chNum += 1;
                }
                while (chNum < ch_filter.Length && ch_filter[chNum] != '@')
                {
                    ch_mz += ch_filter[chNum];
                    chNum += 1;
                }
            }

            double.TryParse(ch_mz, out parent_mz);
            return parent_mz;
        }  // GetParentMz()

        public static int GetMSLevel(IRawDataPlus raw_file, int scan_num)
        {
            var scan_filter = raw_file.GetFilterForScanNumber(scan_num);
            var ms_level = scan_filter.MSOrder;
            if (ms_level < MSOrderType.Ms)  // default to Ms = 1   
            {
                return (int)MSOrderType.Ms;
            }
            return (int)ms_level;
        }

        public static (
            int numPackets,
            double retentionTime,
            double lowMass,
            double highMass,
            double tic,
            double basePeakMass,
            double basePeakIntensity,
            int numChannels,
            bool uniformTime,
            double frequency
            ) GetScanHeaderInfoForScanNum(IRawDataPlus raw_file, int scan_num)
            {
                ScanStatistics scan_statistics = raw_file.GetScanStatsForScanNumber(scan_num);
                return (
                    scan_statistics.PacketCount,
                    scan_statistics.StartTime,
                    scan_statistics.LowMass,
                    scan_statistics.HighMass,
                    scan_statistics.TIC,
                    scan_statistics.BasePeakMass,
                    scan_statistics.BasePeakIntensity,
                    scan_statistics.NumberOfChannels,
                    scan_statistics.IsUniformTime,
                    scan_statistics.Frequency
                );
            } // GetScanHeaderInfoForScanNum()

        public static (double lowMass, double highMass) GetLowHighMassFromFirstScanNumber(IRawDataPlus raw_file, int scan_num)
        {
            ScanStatistics scan_statistics = raw_file.GetScanStatsForScanNumber(scan_num);
            return (
                scan_statistics.LowMass,
                scan_statistics.HighMass
            );
        }  // GetLowHighMassFromFirstScanNumber()

        public static int GetLastSpectrumNumber(IRawDataPlus raw_file)
            {
                return raw_file.RunHeaderEx.LastSpectrum;
            }
            public static int GetNextScanNum(int current_scan_num)
            {
                return(current_scan_num + 1);
            }
            public static double GetScanTime(IRawDataPlus raw_file, int scan_num)
            {
                double start_time = raw_file.RetentionTimeFromScanNumber(scan_num);  // # C++ m_xraw2_class->RTFromScanNum()
                return start_time;
            }

            public static OrderedDictionary GetKeyValuePairs(List<string> keys, List<string> values, int count)
            {
                var result = new OrderedDictionary();
                for (int i = 0; i < count; i++)
                {
                    string key = keys[i].ToString();
                    string value = values[i].ToString();
                    if (result.Contains(key))
                    {
                        string old_val = result[key].ToString();
                        if (!String.IsNullOrEmpty(old_val)) 
                        {
                            // We leave this warning in for our RAW file testing in case we encounter a non-empty value
                            // C++ MprcExtractRaw seems to overwrite the values that its interested in too.
                            Console.WriteLine("WARNING: Removed FROM DICT-- " + key + ": " + old_val);
                        }
                        result.Remove(key);
                    }
                    result.Add(key, value);
                }
                //
                return result;
            }

            public static Tuple<List<string>, List<string>, int> GetTrailerExtraForScanNum(IRawDataPlus rawFile, int scanNum)
            {
                LogEntry trailerInfo = rawFile.GetTrailerExtraInformation(scanNum);
                ScanTrailer scanTrailer = new ScanTrailer(trailerInfo);
                return new Tuple<List<string>, List<string>, int>(scanTrailer.Labels, scanTrailer.Values, scanTrailer.Length);
            }

            public static (double, double, bool, double, double, double, double, double, double, double) GetScanHeaderData(IRawDataPlus rawFile, int scanNum)
            {
                double ionInjectionTimeMs = 0.0;
                double elapsedTimeSeconds = 0.0;
                bool lockMassFound = false;
                double lockMassShift = 0.0;
                double conI = 0.0, conA = 0.0, conB = 0.0, conC = 0.0, conD = 0.0, conE = 0.0;
                Tuple<List<string>, List<string>, int> trailerData = GetTrailerExtraForScanNum(rawFile, scanNum);
                List<string> trailer_data_labels = trailerData.Item1;
                List<string> trailer_data_values = trailerData.Item2;
                int trailer_data_length = trailerData.Item3;
                OrderedDictionary keyValuePairs = GetKeyValuePairs(trailer_data_labels, trailer_data_values, trailer_data_length);
                foreach (string sLabel in keyValuePairs.Keys)
                {
                    string sData = (string)keyValuePairs[sLabel];
                    if (sLabel == "Ion Injection Time (ms):")
                    {
                        ionInjectionTimeMs = double.Parse(sData);
                    }
                    else if (sLabel == "Elapsed Scan Time (sec):")
                    {
                        elapsedTimeSeconds = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter I:")
                    {
                        conI = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter A:")
                    {
                        conA = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter B:")
                    {
                        conB = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter C:")
                    {
                        conC = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter D:")
                    {
                        conD = double.Parse(sData);
                    }
                    else if (sLabel == "Conversion Parameter E:")
                    {
                        conE = double.Parse(sData);
                    }
                    else if (sLabel == "FT Analyzer Message:")
                    {
                        string lockPartStart = "Lock(";
                        string lockPartEnd = ")";
                        string unitSuffix = "ppm";

                        int lockStart = sData.LastIndexOf(lockPartStart);
                        if (lockStart != -1)
                        {
                            int lockEnd = sData.IndexOf(lockPartEnd, lockStart + lockPartStart.Length);
                            if (lockEnd != -1)
                            {
                                string lockInfo = sData.Substring(lockStart + lockPartStart.Length, lockEnd - lockStart - lockPartStart.Length);
                                if (lockInfo.Substring(lockInfo.Length - 3) == unitSuffix)
                                {
                                    if (!lockInfo.Contains("NOT FOUND"))
                                    {
                                        lockMassFound = true;
                                    }
                                }

                                int lmNumStart = lockInfo.LastIndexOf(' ');
                                if (lmNumStart != -1)
                                {
                                    string ppmVal = lockInfo.Substring(lmNumStart + 1, lockInfo.Length - lmNumStart - unitSuffix.Length);
                                    lockMassShift = double.Parse(ppmVal);
                                }
                            }
                        }
                    }
                }
                return (ionInjectionTimeMs, elapsedTimeSeconds, lockMassFound, lockMassShift, conI, conA, conB, conC, conD, conE);
            }  // GetScanHeaderData

            public static (double, List<string>, List<string>, int) GetStatusLogForScanNum(IRawDataPlus rawFile, int scanNum)
            {
                /* Reads and returns the status log data fields present in the RAW file.
                   returns statusLogRT, log_entry.labels, log_entry.values, log_entry.length */
                // Get the status log for this scan
                double statusLogRT = GetScanTime(rawFile, scanNum);
                LogEntry logEntry = rawFile.GetStatusLogForRetentionTime(statusLogRT);
                List<string> logEntryLabels = logEntry.Labels.ToList();
                List<string> logEntryValues = logEntry.Values.ToList();
                int logEntryLength = logEntry.Length;
                return (statusLogRT, logEntryLabels, logEntryValues, logEntryLength);
            }

            public static (double, double, double, double, double, double, double, double, double) GetStatusLogData(IRawDataPlus rawFile, int scanNum)
            {
                /*  This function tries to obtain additional values from the status log
                    :param raw_file:
                    :param scan_num:
                    :return:  sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge, ftVacuumPiraniGauge1,
                    ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp
                */
                double sourceCurrent = 0.0;
                double vacuumIonGauge = 0.0;
                double vacuumConvectronGauge = 0.0;
                double ftVacuumPenningGauge = 0.0;
                double ftVacuumPiraniGauge1 = 0.0;
                double ionMultiplier1 = 0.0;
                double ionMultiplier2 = 0.0;
                double ftCeMeasureVoltage = 0.0;
                double ftAnalyzerTemp = 0.0;
                int nArraySize = 0;
                double statusLogRT = 0.0;
                List<string> logEntryLabels = new List<string>();
                List<string> logEntryValues = new List<string>();
                (statusLogRT, logEntryLabels, logEntryValues, nArraySize) = GetStatusLogForScanNum(rawFile, scanNum);
                OrderedDictionary keyValuePairs = GetKeyValuePairs(logEntryLabels, logEntryValues, nArraySize);
                int section = 0;
                foreach (string sLabel in keyValuePairs.Keys)
                {
                    string sData = (string)keyValuePairs[sLabel];
                    if (sLabel == "API SOURCE" || sLabel == "======  Ion Source:  ======:")
                    {
                        section = 1;
                    }
                    else if (sLabel == "VACUUM" || sLabel == "====== Vacuum: ======:")
                    {
                        section = 2;
                    }
                    else if (sLabel == "FT VACUUM")
                    {
                        section = 3;
                    }
                    else if (sLabel == "TURBO PUMP")
                    {
                        section = 4;
                    }
                    else if (sLabel == "FT TURBO PUMP 1")
                    {
                        section = 5;
                    }
                    else if (sLabel == "FT TURBO PUMP 2")
                    {
                        section = 6;
                    }
                    else if (sLabel == "FT TURBO PUMP 3")
                    {
                        section = 7;
                    }
                    else if (sLabel == "ION OPTICS")
                    {
                        section = 8;
                    }
                    else if (sLabel == "MAIN RF")
                    {
                        section = 9;
                    }
                    else if (sLabel == "ION DETECTION SYSTEM")
                    {
                        section = 10;
                    }
                    else if (sLabel == "FT Analyzer")
                    {
                        section = 11;
                    }
                    else if (sLabel == "POWER SUPPLIES")
                    {
                        section = 12;
                    }
                    else if (sLabel == "FT POWER SUPPLIES")
                    {
                        section = 13;
                    }
                    else if (sLabel == "INSTRUMENT STATUS")
                    {
                        section = 14;
                    }
                    else if (sLabel == "SYRINGE PUMP")
                    {
                        section = 15;
                    }
                    else if (sLabel == "DIVERT VALVE")
                    {
                        section = 16;
                    }

                    // no good switch / case in Python < 3.10
                    if (section == 1) // API SOURCE
                    {
                        if (sLabel == "Source Current (uA):" || sLabel == "Spray Current (µA)")
                        {
                            sourceCurrent = double.Parse(sData);
                        }
                    }
                    else if (section == 2) // VACUUM
                    {
                        if (sLabel == "Ion Gauge (E-5 Torr):")
                        {
                            vacuumIonGauge = double.Parse(sData);
                        }
                        else if (sLabel == "Convectron Gauge (Torr):")
                        {
                            vacuumConvectronGauge = double.Parse(sData);
                        }
                    }
                    else if (section == 3) // FT VACUUM
                    {
                        if (sLabel == "FT Penning Gauge (E-10 Torr):")
                        {
                            ftVacuumPenningGauge = double.Parse(sData);
                        }
                        else if (sLabel == "FT Pirani Gauge 1 (Torr):")
                        {
                            ftVacuumPiraniGauge1 = double.Parse(sData);
                        }
                    }
                    else if (section == 10) // ION DETECTION SYSTEM
                    {
                        if (sLabel == "Multiplier 1 (V):")
                        {
                            ionMultiplier1 = double.Parse(sData);
                        }
                        else if (sLabel == "Multiplier 2 (V):")
                        {
                            ionMultiplier2 = double.Parse(sData);
                        }
                    }
                    else if (section == 11) // FT Analyzer
                    {
                        string sLabelCompared = sLabel.Substring(0, "FT Analyzer Temp. (".Length);  // check this 
                        string sLabelComparing = "FT Analyzer Temp. (";
                        if (sLabel == "FT CE Measure Voltage (V):")
                        {
                            ftCeMeasureVoltage = double.Parse(sData);
                        }
                        else if (sLabelCompared == sLabelComparing)
                        {
                            ftAnalyzerTemp = double.Parse(sData);
                        }
                    }
                }  // foeach
                return (sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge, ftVacuumPiraniGauge1,
                        ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp);
        }  //  GetStatusLogData

        public static (double, double, double, double, int, double, double, double, 
            double, bool, double, double, double, double, double, double, double, 
            double, double, double, double, double, 
            double, double, double, double, double, double)
        /*
            Obtain information about the scan
            :param raw_file:
            :param scan_num:
            :return:   tic, retentionTime, lowMass, highMass, childScans, ionInjectionTimeMs, cycleTimeSeconds elapsedTimeSeconds,
            timeToNextScanSeconds, lockMassFound, lockMassShift, conI, conA, conB, conC, conD, conE,
            sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge, ftVacuumPiraniGauge1,
            ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp, basePeakMass, basePeakIntensit
         */
        // numPackets:  How many mass/intensity pairs
        // numChannels: Number of channels acquired for the scan
        // uniformTime: Indicating whether the sampling time increment for the current controller is uniform (?)
        // frequency:   Sampling frequency for the current controller
        GetScanInfo(IRawDataPlus raw_file, int scan_num)
        {
            int numPackets;
            double retentionTime;
            double lowMass;
            double highMass;
            double tic;
            double basePeakMass;
            double basePeakIntensity;
            int numChannels;
            bool uniformTime;
            double frequency;
            (numPackets, retentionTime, lowMass, highMass, tic, basePeakMass, basePeakIntensity, numChannels, uniformTime, frequency) =
                GetScanHeaderInfoForScanNum(raw_file, scan_num);
            // Extract delta information by parsing subsequent scans
            double timeToNextScanSeconds = 0.0;
            double cycleTimeSeconds = 0.0;
            int elapsedTimeSeconds = 0;
            int childScans = 0;
            if (GetMSLevel(raw_file, scan_num) == 1)
            {
                int nextScan = GetNextScanNum(scan_num);
                while (nextScan <= GetLastSpectrumNumber(raw_file))
                {
                    int msLevel = GetMSLevel(raw_file, nextScan);
                    double nextRT = GetScanTime(raw_file, nextScan);
                    if (msLevel == 1)
                    {
                        cycleTimeSeconds = (nextRT - retentionTime) * 60.0;
                        break;
                    }
                    else
                    {
                        childScans++;
                    }
                    nextScan = GetNextScanNum(nextScan);
                }
            }

            if ((scan_num + 1) <= GetLastSpectrumNumber(raw_file))
            {
                double nextRT = GetScanTime(raw_file, scan_num + 1);
                timeToNextScanSeconds = (nextRT - retentionTime) * 60.0;
            }

            double ionInjectionTimeMs;
            double elapsedTimeSecondsHeader;
            bool lockMassFound;
            double lockMassShift;
            double conI;
            double conA;
            double conB;
            double conC;
            double conD;
            double conE;
            (ionInjectionTimeMs, elapsedTimeSecondsHeader, lockMassFound, lockMassShift,
                conI, conA, conB, conC, conD, conE) = GetScanHeaderData(raw_file, scan_num);

            double sourceCurrent;
            double vacuumIonGauge;
            double vacuumConvectronGauge;
            double ftVacuumPenningGauge;
            double ftVacuumPiraniGauge1;
            double ionMultiplier1;
            double ionMultiplier2;
            double ftCeMeasureVoltage;
            double ftAnalyzerTemp;
            (sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge,
                ftVacuumPiraniGauge1, ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp) = 
                GetStatusLogData(raw_file, scan_num); 

            return (tic, retentionTime, lowMass, highMass, childScans, ionInjectionTimeMs, cycleTimeSeconds, elapsedTimeSeconds,
                     timeToNextScanSeconds, lockMassFound, lockMassShift, conI, conA, conB, conC, conD, conE,
                     sourceCurrent, vacuumIonGauge, vacuumConvectronGauge, ftVacuumPenningGauge, ftVacuumPiraniGauge1,
                     ionMultiplier1, ionMultiplier2, ftCeMeasureVoltage, ftAnalyzerTemp, basePeakMass, basePeakIntensity);
        }  // GetScanInfo()

        public static int GetParentScan(IRawDataPlus raw_file, int scan_num)
        {
            int msN_level = GetMSLevel(raw_file, scan_num);
            int level = msN_level;
            int i = 0;
            while (level >= msN_level && (scan_num - i) >= 1)
            {
                level = GetMSLevel(raw_file, scan_num - i);
                i++;
            }
            i--;
            return scan_num - i;
        }

        // # Fills the buf with dissociation type. The buffer must be at least 10 characters long (+trailing zero)
        public static string GetDissociationType(IRawDataPlus raw_file, int scanNum)
        {
            // Many ways to skin this cat from the C++ GetFilterForScanNum(), scan_statistics.scan_type seems simplest
            ScanStatistics scanStatistics = raw_file.GetScanStatsForScanNumber(scanNum);
            string scanFilterName = scanStatistics.ScanType;
            int inIdx = scanFilterName.IndexOf('@');
            if (inIdx == -1)
            {
                return string.Empty;
            }
            inIdx += 1;
            int outIdx = 0;
            string output = string.Empty;
            while (inIdx < scanFilterName.Length && !char.IsDigit(scanFilterName[inIdx]) && outIdx < 10)
            {
                output += scanFilterName[inIdx];
                inIdx += 1;
                outIdx += 1;
            }
            return output;
        }  // GetDissociationType()

        public static (double,double) GetBasePeak(double[] mzs, double[] intensities, double basePeakMz, double basePeakIntensity, double minDistanceSecondFromBaseDa)
        {
            double secondPeakMz = 0.0;
            double secondPeakIntensity = 0.0;
            for (int i = 0; i < mzs.Length; i++)
            {
                double mzsIter = mzs[i];
                double intIter = intensities[i];
                if (intIter > secondPeakIntensity && Math.Abs(mzsIter - basePeakMz) > minDistanceSecondFromBaseDa)
                {
                    secondPeakIntensity = intIter;
                    secondPeakMz = mzsIter;
                }
            }
            return (secondPeakMz, secondPeakIntensity);
        }  // GetBasePeak()

        public static void extractInstrumentMethodData(IRawDataPlus rawFile, string instrumentMethodFileName)
        {
            MetadataWriter metadataWriter = new MetadataWriter(instrumentMethodFileName);
            metadataWriter.WriteMetadata(rawFile);
            Console.WriteLine("INFO: Finished parsing for instrument method data " + instrumentMethodFileName);
        }

    }  // class mPRCUtils
}  // namespace

