/*
    MprcExtractRaw Globals live in here.
*/
// Version of MprcExtractRaw3
namespace MprcExtractRaw3
{
    using System.Collections.Generic;
    public static class Config
    {
        // Version of MprcExtractRaw3
        public const string MPRC_EXTRACTRAW3_VERSION = "3.00.00";

        // Polymer detection
        public const int MIN_SEGMENT_SIZE = 30;
        public const int MAX_SEGMENT_SIZE = 100;

        // Command line options
        public const string DATA = "--data";
        public const string MZ_RANGE = "--mzrange";
        public const string PARAMS_FILE = "--params";
        // MPRC params
        public const string RAW_FILE = "--raw";
        // Data params
        public const string INFO_FILE = "--info";
        public const string SPECTRA_FILE = "--spectra";
        public const string CHROMATOGRAM_FILE = "--chromatogram";
        public const string TUNE_METHOD_FILE = "--tune";
        public const string INSTRUMENT_METHOD_FILE = "--instrument";
        public const string SAMPLE_INFORMATION_FILE = "--sample";
        public const string ERROR_LOG_FILE = "--errorlog";
        public const string UV_DATA_FILE = "--uv";  // From UV controller (pressure info)
        public const string ANALOG_DATA_FILE = "--analog";  // From Analog / Digital controller (pressure info for new instruments)
        public const string RTC_FILE = "--rtc";  // File to extract retention time calibration data (e.g. Pierce)
        public const string RTC_PRECURSOR_MZS = "--rtcPrecursorMzs";  // Colon-separated list of precursor m/z values
        public const string RTC_MASS_ROUNDING = "--rtcMassRounding";  // number of decimals for rounding the m/z values, 2 by default
        public const string RTC_PPM_MASS_TOL = "--rtcPpmMassTol";  // mass tolerance for looking up the precursor m/z values in the spectral data
                                                                   // MZ range params (rawFile is shared with MPRC params)
        public const string MIN_MZ = "--min";
        public const string MAX_MZ = "--max";
        public const string PEAKS_FILE = "--peaks";
        public const string FIND_POLYMERS_NUM_PROCESSORS = "--findPolymersNumProcessors";
        public const string DEBUG = "--debug";

        // Settings
        public const int CHROMATOGRAM_MZ_BINS = 1000;  // How many bins per m/z range to report
        public const double SECOND_PEAK_MIN_DISTANCE_FROM_BASE = 5.0;  // How far can the second most abundant peak be from the base peak (Da)

        // Column names for the spectra data file
        public const string SCAN_ID = "Scan Id";
        public const string PARENT_MZ = "Parent m/z";
        public const string TOTAL_ION_CURRENT = "TIC";
        public const string RETENTION_TIME = "RT";
        public const string MS_LEVEL = "MS Level";
        public const string PARENT_SCAN = "Parent Scan";
        public const string CHILD_SCANS = "Child Scans";
        public const string ION_INJECTION_TIME_MS = "Ion Injection Time";
        public const string CYCLE_TIME_SECONDS = "Cycle Time";
        public const string ELAPSED_TIME_SECONDS = "Elapsed Time";
        public const string DEAD_TIME_SECONDS = "Dead Time";
        public const string TIME_TO_NEXT_SCAN_SECONDS = "Time To Next Scan";
        public const string LOCK_MASS_FOUND = "Lock Mass Found";
        public const string LOCK_MASS_SHIFT = "Lock Mass Shift";
        public const string CON_I = "Conversion Parameter I";
        public const string CON_A = "Conversion Parameter A";
        public const string CON_B = "Conversion Parameter B";
        public const string CON_C = "Conversion Parameter C";
        public const string CON_D = "Conversion Parameter D";
        public const string CON_E = "Conversion Parameter E";
        public const string DISSOCIATION_TYPE = "Dissociation Type";
        public const string POLYMER_SEGMENT = "Polymer Segment Size";
        public const string POLYMER_OFFSET = "Polymer Offset";
        public const string POLYMER_SCORE = "Polymer Score";
        public const string POLYMER_PVALUE = "Polymer p-value";
        public const string POLYMER_FOR_MASS_OFFSET = "Polymer ({0}) Offset";
        public const string POLYMER_FOR_MASS_SCORE = "Polymer ({0}) Score";
        public const string POLYMER_FOR_MASS_PVALUE = "Polymer ({0}) p-value";
        public const string BASE_PEAK_MZ = "Base Peak m/z";
        public const string BASE_PEAK_INTENSITY = "Base Peak Intensity";
        public const string SECOND_PEAK_MZ = "Second Peak m/z";
        public const string SECOND_PEAK_INTENSITY = "Second Peak Intensity";

        // API SOURCE
        public const string SOURCE_CURRENT = "Source Current (uA)";
        public const string SPRAY_CURRENT = "Spray Current (uA)";
        // VACUUM
        public const string VACUUM_ION_GAUGE = "Vacuum Ion Gauge (E-5 Torr)";
        public const string VACUUM_CONVECTRON_GAUGE = "Vacuum Convectron Gauge (Torr)";
        // FT VACUUM
        public const string FT_VACUUM_PENNING_GAUGE = "FT Penning Gauge (E-10 Torr)";
        public const string FT_VACUUM_PIRANI_GAUGE1 = "FT Pirani Gauge 1 (Torr)";
        // ION DETECTION SYSTEM
        public const string ION_MULTIPLIER1 = "Multiplier 1 (V)";
        public const string ION_MULTIPLIER2 = "Multiplier 2 (V)";
        // FT ANALYZER
        public const string FT_CE_MEASURE_VOLTAGE = "FT CE Measure Voltage (V)";
        public const string FT_ANALYZER_TEMP = "FT Analyzer Temp (C)";

        // Parameter names for raw info data file
        public const string ORIGINAL_RAW_FILE_NAME = "Name";
        public const string CREATION_DATE = "Creation Date";
        public const string SAMPLE_ID = "Sample Id";
        public const string NUMBER_OF_MS1_SPECTRA = "MS1 Spectra";
        public const string NUMBER_OF_MS2_SPECTRA = "MS2 Spectra";
        public const string NUMBER_OF_MS3PLUS_SPECTRA = "MS3+ Spectra";
        public const string INSTRUMENT_NAME = "Instrument Name";
        public const string INSTRUMENT_SERIAL = "Instrument Serial";
        public const string RUN_TIME_IN_SECONDS = "Run Time (seconds)";
        public const string COMMENT = "Comment";

        // List of extra polymer masses to investigate
        public static readonly List<int> POLYMER_MASSES = new List<int> { 162 };
    }

}
