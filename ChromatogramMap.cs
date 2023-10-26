namespace MprcExtractRaw3
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;
    using System.Threading.Tasks;
    using System.Drawing;
    using System.Drawing.Imaging;
    using System.Runtime.Serialization;
    using System.IO;
    using System.Text.RegularExpressions;
    using System.Runtime.InteropServices;

    // re-enable when magick.NET does not crash on odin (Centos7)
    //using ImageMagick;
    using System.Xml.Linq;

    /// <summary>
    /// Class used for direct memory access to 8bit grayscale images.
    /// Not thread-safe.
    /// ref https://www.sharpgis.net/post/Working-with-8bit-images-in-NET
    /// </summary>
    public class Image8Bit : IDisposable
    {
        private BitmapData bmd;
        private Bitmap b;
        /// <summary>
        /// Locks an 8bit image in memory for fast get/set pixel functions.
        /// Remember to Dispose object to release memory.
        /// 
        /// </summary>
        /// Bitmap reference
        public Image8Bit(Bitmap bitmap)
        {
            if (bitmap.PixelFormat != System.Drawing.Imaging.PixelFormat.Format8bppIndexed)
                throw (new System.Exception("Invalid PixelFormat. 8 bit indexed required"));
            b = bitmap; //Store a private reference to the bitmap
            bmd = b.LockBits(new Rectangle(0, 0, b.Width, b.Height),
                                ImageLockMode.ReadWrite, b.PixelFormat);
        }

        /// <summary>
        /// Releases memory
        /// </summary>
        public void Dispose()
        {
            b.UnlockBits(bmd);
        }

        /// <summary>
        /// Gets color of an 8bit-pixel
        /// </summary>
        /// <param name="x">Row</param>
        /// <param name="y">Column</param>
        /// <returns>Color of pixel</returns>
        public unsafe System.Drawing.Color GetPixel(int x, int y)
        {
            byte* p = (byte*)bmd.Scan0.ToPointer();
            //always assumes 8 bit per pixels
            int offset = y * bmd.Stride + x;
            return GetColorFromIndex(p[offset]);
        }

        /// <summary>
        /// Sets color of an 8bit-pixel
        /// </summary>
        /// <param name="x">Row</param>
        /// <param name="y">Column</param>
        /// <param name="c">Color index</param>
        public unsafe void SetPixel(int x, int y, byte c)
        {
            byte* p = (byte*)bmd.Scan0.ToPointer();
            //always assumes 8 bit per pixels
            int offset = y * bmd.Stride + (x);
            p[offset] = c;
        }

        /// <summary>
        /// Sets the palette for the referenced image to Grayscale
        /// </summary>
        public void MakeGrayscale()
        {
            SetGrayscalePalette(this.b);
        }

        /// <summary>
        /// Sets the palette of an image to grayscales (0=black, 255=white)
        /// </summary>
        /// <param name="b">Bitmap to set palette on</param>
        public static void SetGrayscalePalette(Bitmap b)
        {
            ColorPalette pal = b.Palette;
            for (int i = 0; i < 256; i++)
                pal.Entries[i] = Color.FromArgb(255, i, i, i);
            b.Palette = pal;
        }

        private System.Drawing.Color GetColorFromIndex(byte c)
        {
            return b.Palette.Entries[c];
        }
    }

    public class ChromatogramMap
    {
        private const int MIN_INTENSITY = 1;
        private const int DETAILED_HISTOGRAM_BUCKETS = 65536;
        private const int EXPORT_HISTOGRAM_BUCKETS = 256;

        private int massResolution;
        private double minMass;
        private double maxMass;
        private double maxLogIntensity;
        private double[] data;
        private int dataNumRows;

        public ChromatogramMap(int massResolution, double minMass, double maxMass)
        {
            this.massResolution = massResolution;  // How many pixels per the width of a spectrum
            this.minMass = minMass;  // Minimum mass to capture
            this.maxMass = maxMass;  // Maximum mass to capture
            this.maxLogIntensity = 0.0;  // float maximum log (peak intensity)
            // Condensed data of all spectra so far. log of each value is being stored. Values <1 are stored as 0
            // internal rep is a 1D list indexed by row*num_cols + col, massResolution is the length of a row
            this.data = new double[0];
            this.dataNumRows = 0;
        }

        ~ChromatogramMap()
        {
            this.data = null;
        }

        public void SetData(double[] values)
        {
            this.data = this.data.Concat(values).ToArray();
            this.dataNumRows++;
        }

        // SetData multithreaded version sets data en masse
        public void SetData(double[] values, int num_rows)
        {
            this.data = this.data.Concat(values).ToArray();
            this.dataNumRows += num_rows;
        }

        public double[] GetData()
        {
            return this.data;
        }

        public void SetMaxLogIntensity(double maxLogIntensity)
        {
            this.maxLogIntensity = maxLogIntensity;
        }
        public double GetMaxLogIntensity()
        {
            return this.maxLogIntensity;
        }

        // Add spectrum to the in-memory representation
        public void AddSpectrum(double[] inMzs, double[] inIntensities)
        {
            double[] values = new double[this.massResolution];
            double MASSRANGE = this.maxMass - this.minMass;
            int mzsLen = inMzs.Length;
            double[] minMasses = new double[mzsLen];
            double[] MASSRANGES = new double[mzsLen];
            double[] massResolutions = new double[mzsLen];
            for (int idx = 0; idx < mzsLen; idx++)
            {
                minMasses[idx] = this.minMass;
                MASSRANGES[idx] = MASSRANGE;
                massResolutions[idx] = this.massResolution;
            }
            //  Where does the value fit into our array
            double[] buckets = new double[mzsLen];
            int[] bucketIndices = new int[mzsLen];
            for (int i = 0; i < mzsLen; i++)
            {
                buckets[i] = ((inMzs[i] - minMasses[i]) / MASSRANGES[i]) * massResolutions[i];
                bucketIndices[i] = (int)buckets[i];
            }
            // Intensity is log-scaled, has to be >= 1
            double[] intensities = new double[inIntensities.Length];
            for (int i = 0; i < inIntensities.Length; i++)
            {
                intensities[i] = inIntensities[i] < MIN_INTENSITY ? MIN_INTENSITY : inIntensities[i];
            }
            double[] logIntensities = new double[intensities.Length];
            for (int i = 0; i < intensities.Length; i++)
            {
                logIntensities[i] = Math.Log(intensities[i]);
            }
            for (int i = 0; i < logIntensities.Length; i++)
            {
                double logIntensity = logIntensities[i];
                int bucket = bucketIndices[i];
                if (bucket >= 0 && bucket < this.massResolution)
                {
                    if (logIntensity > this.maxLogIntensity)
                        this.maxLogIntensity = logIntensity;
                    // Intensity is log-scaled, has to be >= 1
                    if (values[bucket] < logIntensity)
                        values[bucket] = logIntensity;
                }
            }
            SetData(values);
        }

        // Output the equalized picture into a given file in .gif or .png format
        public void DumpEqualized(string filename)
        {
            string ext = Path.GetExtension(filename);
            string output_dir = Directory.GetParent(filename).FullName + Path.DirectorySeparatorChar;
            string base_filename = Path.GetFileNameWithoutExtension(filename);
            // error check if we have no data i.e. spectrum
            if (GetHeight() == 0)
            {
                // output yellow warning text with black background
                Bitmap errorImage = new Bitmap(256, 256, PixelFormat.Format24bppRgb);
                Graphics graphics = Graphics.FromImage(errorImage);
                graphics.FillRectangle(new SolidBrush(Color.Black), 0, 0, 256, 256);
                Font font = SystemFonts.DefaultFont;
                graphics.DrawString("No spectrum was found", font, new SolidBrush(Color.Yellow), 10, 10);
                //string er 
                ImageFormat imgFormat = ImageFormat.Gif;
                errorImage.Save(filename, imgFormat);
                graphics.Dispose();
                errorImage.Dispose();
                // BEG TMP CODE
                // string errorComment = "mz:" + this.minMass.ToString() + "," + this.maxMass.ToString();
                //using (MagickImage image2 = new MagickImage(filename))
                //{
                //    image2.Comment = errorComment;
                //    image2.Write(filename);
                //}
                // END TMP CODE
                return;
            }

            int[] histogramData = new int[DETAILED_HISTOGRAM_BUCKETS];
            Histogram(histogramData, DETAILED_HISTOGRAM_BUCKETS);
            double[] exportHistogram = new double[EXPORT_HISTOGRAM_BUCKETS + 1];
            Equalize(histogramData, DETAILED_HISTOGRAM_BUCKETS, GetWidth() * GetHeight(), 
                this.maxLogIntensity, exportHistogram, EXPORT_HISTOGRAM_BUCKETS);

            int outputWidth = GetHeight();
            int outputHeight = GetWidth();
            Bitmap image = new Bitmap(outputWidth, outputHeight, PixelFormat.Format8bppIndexed);
            Bitmap image_I24B = new Bitmap(outputWidth, outputHeight, PixelFormat.Format24bppRgb);
            // Set up the color palette
            ColorPalette palette = image.Palette;
            int[] colorMap = new int[3 * (EXPORT_HISTOGRAM_BUCKETS + 1)];
            palette.Entries[0] = Color.FromArgb(0, 0, 0);
            for (int i = 1; i < EXPORT_HISTOGRAM_BUCKETS; i++)
            {
                // Histogram is log-scale.
                // Convert to actual value
                double bucketCentroid = (Math.Exp(exportHistogram[i - 1]) + Math.Exp(exportHistogram[i])) / 2.0;
                // Convert to log-10 scale - 1 is zero, 1E10 is 1
                double log10scale = Math.Log10(bucketCentroid) / 10.0;
                log10scale = Math.Max(0.0, Math.Min(1.0, log10scale));
                // Convert to 3-byte value
                long color = (long)(0xFFFFFF * log10scale);
                // R,G,B form together intensity, log scale, 255, 255, 255 being the maximum intensity (1E10)
                byte red = (byte)((color & 0xFF0000) >> 16);
                byte green = (byte)((color & 0x00FF00) >> 8);
                byte blue = (byte)(color & 0x0000FF);
                colorMap[3 * i] = red;
                colorMap[3 * i + 1] = green;
                colorMap[3 * i + 2] = blue;
                palette.Entries[i] = Color.FromArgb(red, green, blue);  // used for 8-bit images
            }
            image.Palette = palette;
            Image8Bit image_I8B = new Image8Bit(image);

            // Now we can go value by value, get the resulting buckets and output a byte
            // Title: Inj19_PooledAmyloid10_20220824_Renal_Optimizedv1_S1.raw.chroma.gif
            // assuming orig order was left => right and top => bottom
            int outputRow = 0;
            for (int values = GetWidth() - 1; values >= 0; values--)
            {
                int outputCol = 0;
                for (int spectrum = 0; spectrum < GetHeight(); spectrum++)
                {
                    double intensity = this.data[spectrum * GetWidth() + values];   
                    int bucket = GetBucket(intensity, exportHistogram, EXPORT_HISTOGRAM_BUCKETS);
                    int colorMapIdx = Math.Max(0, Math.Min(255, bucket));
                    // pixels[output_row][output_col] = np.ubyte(bucket)  # np indexed by row, col
                    image_I8B.SetPixel(outputCol, outputRow, (byte)colorMapIdx);

                    colorMapIdx *= 3;  //  index into the start of our red, green, blue in the colormap
                    Color color = Color.FromArgb(colorMap[colorMapIdx], colorMap[colorMapIdx + 1], colorMap[colorMapIdx + 2]);
                    // pixels[output_row][output_col] = np.ubyte(bucket)  # np indexed by row, col
                    image_I24B.SetPixel(outputCol, outputRow, color);
                    outputCol++;
                }
                outputRow++;
            }
            // We save the 8-bit GIF for R's caTools in order to extract its palette.
            // We save the 24-bit GIF for R's magick, because 8-bit GIFs appear blank under Mono (Centos7).
            // BUG blnak GIFs mono-- See https://github.com/mono/mono/issues/12425

            image.Save(filename, ImageFormat.Gif);
            string filename_I24B = output_dir + base_filename + "_I24B"  + ext;
            image.Save(filename_I24B, ImageFormat.Gif);

            // Add comment to image file as it is used in swift/scripts/src/main/resources/bin/util/rPpmPlot.r
            // ref https://github.com/python-pillow/Pillow/issues/6024 maybe try translating to C# via G and try if still needed
            string comment = "mz:" + this.minMass.ToString() + "," + this.maxMass.ToString();
            
            // re-enable when magick.NET does not crash on odin (Centos7)
            //using (MagickImage image2 = new MagickImage(filename))
            //{
            //    image2.Comment = comment;
            //    image2.Write(filename);
            //}
       
            // "CR-23-28069_Delm_20230520_QE2_s1.raw.chroma_mzdims.tsv"
            string mzdims_filename = output_dir + base_filename + "_mzdims" + ".tsv";
            try
            {
                using (StreamWriter writer = new StreamWriter(mzdims_filename))
                {
                    writer.WriteLine("minMz\tmaxMz");
                    writer.WriteLine($"{this.minMass}\t{this.maxMass}");
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("ERROR: Unable to create chromatogram M/Z dimension information in DumpEqualized(): " + ex.Message);
                Environment.Exit(1);
            }

            image_I8B.Dispose();
            image.Dispose();
            image_I24B.Dispose();
        }

        public static bool DumpEqualizedWrapper(string output_dir, string chroma_base_filename, string chroma_output_filename,
           string FILE_EXT, int massResolution, double lowMass, double highMass, double maxLogIntensity, bool debug)
        {
            ChromatogramMap map = new ChromatogramMap(massResolution, lowMass, highMass);
            map.SetMaxLogIntensity(maxLogIntensity);

            string[] chroma_files = Directory.GetFiles(output_dir, chroma_base_filename + ".*" + FILE_EXT);
            // Sort the filenames based on the extracted number part ex. 200ngHela_20230515_Ex1_1sAuto_SortHighestMz_Inj27.raw.chroma.gif
            string pattern = @"(?<=raw\.chroma\.)\d+(?=\.tsv)";
            var sorted_chroma_files = chroma_files.OrderBy(filename =>
            {
                Match match = Regex.Match(filename, pattern);
                if (match.Success && int.TryParse(match.Value, out int number))
                {
                    return number;
                }
                return int.MaxValue; // Invalid numbers are treated as maximum value
            });

            if (sorted_chroma_files.Count() == 0)
            {
                Console.WriteLine("ERROR: No Chromatogram Map files were generated for {0}", chroma_output_filename);
                return false;
            }
            foreach (string file in sorted_chroma_files)
            { 
                // Using a list may be slow?
                var dataList = new System.Collections.Generic.List<double>();
                using (StreamReader reader = new StreamReader(file))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        string[] elements = line.Split('\t');
                        foreach (string element in elements)
                        {
                            double value = double.Parse(element);
                            dataList.Add(value);
                        }
                    }  // while read each line on each sorted chroma file
                }
                // Convert the list to a double array with the correct size
                double[] dataArray = dataList.ToArray();
                int num_rows = dataArray.Length / massResolution;
                map.SetData(dataArray, num_rows);
            }  // foreach sorteed chroma file

            map.DumpEqualized(chroma_output_filename);

            // cleanup
            if (!debug)
            {
                foreach (string file in sorted_chroma_files)
                {
                    try
                    {
                        System.IO.File.Delete(file);
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("ERROR: deleting the temporary Chroma file: {0}", file);
                        Console.WriteLine(ex.Message);
                    }
                }
            }  // if not debug

            return true;
        }


        public int GetHeight()
        {
            // The number of spectra.
            return this.dataNumRows;
        }

        public int GetWidth()
        {
            //  The number of pixels per spectrum or mass resolution.
            return this.massResolution;
        }

        /*
            Calculate a histogram of intensity values of given size.
            The interval <0, maxIntensity> is split into numBuckets equal buckets.
            We return how many values fall into a given bucket
            This function has been vectorized via numpy
        */ 
        private void Histogram(int[] buckets, int numBuckets)
        {
            int width = GetWidth();
            int height = GetHeight();
            for (int spectrum = 0; spectrum < height; spectrum++)
            {
                int startIdx = spectrum * width;
                int endIdx = startIdx + width;
                double[] intensities = new double[width];
                // intensities = self.data[start_idx:end_idx]
                Array.Copy(this.data, startIdx, intensities, 0, width);
                double[] maxLogIntensities = new double[width];
                int[] numBucketsArray = new int[width];
                for (int i = 0; i < width; i++)
                {
                    maxLogIntensities[i] = this.maxLogIntensity;
                    numBucketsArray[i] = numBuckets;
                }
                double[] bucketIndices = new double[width];
                for (int i = 0; i < width; i++)
                {
                    bucketIndices[i] = (intensities[i] / maxLogIntensities[i]) * numBucketsArray[i];
                }
                int[] bucketIndicesInt = new int[width];
                for (int i = 0; i < width; i++)
                {
                    bucketIndicesInt[i] = (int)bucketIndices[i];
                    if (bucketIndicesInt[i] >= numBuckets)
                    {
                        bucketIndicesInt[i] = numBuckets - 1;
                    }
                }
                for (int i = 0; i < width; i++)
                {
                    buckets[bucketIndicesInt[i]]++;
                }
            }
        }  // Histogram()

        private static void Equalize(int[] inputBuckets, int numInputBuckets, int totalValues, double maxInputIntensity, 
            double[] outputBucketBoundaries, int numOutputBuckets)
        {
            /*
                Equalize a calculated histogram for given output amount of buckets
                :param numInputBuckets: numInputBuckets - original histogram
                :param totalValues: total amount of entries in the original histogram (could be obtained by summing the entire histogram).
                :param maxInputIntensity: maximum intensity of the input histogram (minimum is known to == 0)
                :param outputBucketBoundaries: boundaries of the output histogram buckets (preallocated). There is one extra boundary at the end of the last bucket.
                :param numOutputBuckets: how many output buckets to produce
                :return:
            */
            int consumedValues = 0;
            outputBucketBoundaries[0] = 0.0;
            int inputBucket = 0;
            for (int i = 0; i < numOutputBuckets; i++)
            {
                int fill = 0;
                int expectedFill = (totalValues - consumedValues) / (numOutputBuckets - i);
                while (fill < expectedFill && inputBucket < numInputBuckets)
                {
                    int current = inputBuckets[inputBucket];
                    consumedValues += current;
                    fill += current;
                    inputBucket++;
                }
                double lastInputBucketBoundary = maxInputIntensity * inputBucket / numInputBuckets;
                outputBucketBoundaries[i + 1] = lastInputBucketBoundary;
            }
        }

        private int GetBucket(double input, double[] bucketBoundaries, int numBuckets)
        {
            /*
                Binary histogram search, returns the histogram bucket containing the given value.
                bucket[n] = <bucketBoundaries[n], bucketBoundaries[n+1])
                :param bucketBoundaries:  define buckets (see above)
                :param numBuckets:   is one less than number of bucket boundaries
                :return: the bucket number encoded
            */
            int left = 0;
            int right = numBuckets;
            while (right - left > 1)
            {
                int mid = (left + right) / 2;
                if (bucketBoundaries[mid] >= input)
                {
                    right = mid;
                }
                else
                {
                    left = mid;
                }
            }
            return left;
        }  // GetBucket()

        public void SavePartialChromaMapFile(string filename)
        {
            // internal rep is a 1D list indexed by row*num_cols + col, massResolution is the length of a row
            double[] data = GetData();
            using (StreamWriter writer = new StreamWriter(filename))
            {
                for (int i = 0; i < data.Length; i++)
                {
                    // Convert the double value to string with 12 bits of precision
                    string formattedValue = data[i].ToString("F12");
                    writer.Write(formattedValue);
                    // Add a tab character after each element (except for the last one)
                    if (i < data.Length - 1)
                    {
                        writer.Write("\t");
                    }
                }
            }  // using stream writer
        }  // SavePartialChromaMapFile

    } // calss ChromatogramMap

}
