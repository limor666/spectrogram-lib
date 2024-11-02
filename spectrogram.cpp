#include <cstdio>
#include <cstdint>
#include <cmath> // log()
#include <cfloat> // FLT_MAX, FLT_MIN
#include <cstring> // memset()

#include "colormaps.h" // plasma and such
#include "filterwindow.cpp" // Hanning, Hamming..
#include "mel.cpp" // log-distort image

// rdft uses double but one could #define double as float and compile as such.
// here we just copy to new array of doubles since we must anyway multiply audio samples by array of filter Window.
// Ooura library doesn't have an include file with function definitions
void rdft(int, int, double *, int *, double *); // N, direction, inout, private, sin\cos


class spectrogram
{
    public:
    float *audio_data = nullptr; // pcm = pulse code modulation
    int audio_samples = 0;
    float audio_min = FLT_MAX; // lowest sample value
    float audio_max = FLT_MIN; // highest

    float *image_data = nullptr; // spectrogram data
    int image_lines = 0; // number of fft's
    int image_width = 0; // number of frequencies per fft (half length of sound sample window)

    float *mel_image_data = nullptr; // mel's distortion of spectrogram image 
    int mel_image_width = 0; // number of fft frequencies "bins"

    float *rgb_image_data = nullptr; // colored fft or mel image;
    bool rgb_image_fft_or_mel;

    ~spectrogram() {
        // audio_data buffer is created outside this class
        if (image_data) delete[] image_data;
        if (mel_image_data) delete[] mel_image_data;
        if (rgb_image_data) delete[] rgb_image_data;
    }
    spectrogram(float *audio_samples_float32, int audio_nsamples)
    {
        audio_data = audio_samples_float32;
        audio_samples = audio_nsamples;

        float *p = audio_data;
        for (int i=0; i<audio_nsamples; i++)
        {
            float x = *p++;
            if (x>audio_max) audio_max = x;
            else if (x<audio_min) audio_min = x;
        }
    }

    void normalize_audio(float min=-1, float max=1) {
        float *p = audio_data;
        float diff1 = audio_max-audio_min;
        float diff2 = max-min;
        for (int i=0; i<audio_samples; i++) {
            *p = diff2*(*p - audio_min)/diff1+min;
            p++;
        }
        audio_min=min;
        audio_max=max;
    }

    void spectrogram_from_audio(int start, int stride, int fftsize, fft_window::Type window_type) {
        int WINDOW = fftsize;
        image_width = WINDOW/2; // save in class data
        int NMAX = WINDOW;
        const int NMAXSQRT = 64;
        const int SAMPLESIZE = sizeof(float);

        // buffers prerequired by Ooura FFT
        int ooura_ip[NMAXSQRT + 2];
        ooura_ip[0] = ooura_ip[1] = 0; // first time only for rdft
        double ooura_w[NMAX];
        double ooura_fft_buf[NMAX + 1];

        fft_window window = fft_window(WINDOW, window_type); // create Hanning filter
 
        // first_line window must start at sample > -window
        // last_line window must start at sample <= last_sample
        if (start <= -WINDOW) start = -WINDOW+1;
        image_lines = (audio_samples-start+stride-1)/stride;
        image_data = new float[image_lines*image_width];
    
        double maxfft=-DBL_MAX;
        double minfft=DBL_MAX;
 

        float *img = image_data;
        for (int l=0; l<image_lines; l++) {

            // Multiply audio window by Hamming/Hanning window
            // Consider cases where window may start before audio or end after audio
            int i = l*stride + start;
            if (i < 0) {
                int wl = WINDOW + i;
                memset(ooura_fft_buf,0,-i);
                double *f = ooura_fft_buf-i;
                float *a = audio_data;
                double *w = window.data-i;
                for (int j=0; j < wl; j++) *f++ = *a++ * *w++;
            } else if (i > audio_samples-WINDOW) {
                int wl = audio_samples-i;
                memset(ooura_fft_buf+wl,0, WINDOW-wl);
                double *f = ooura_fft_buf;
                float *a = audio_data+i;
                double *w = window.data;
                for (int j=0; j < wl; j++) *f++ = *a++ * *w++;               
            } else {
                int wl = WINDOW;
                //memset(ooura_fft_buf+wl,0, WINDOW-wl);
                double *f = ooura_fft_buf;
                float *a = audio_data+i;
                double *w = window.data;
                for (int j=0; j < wl; j++) *f++ = *a++ * *w++;
            }

            // Run the Ooura FFT            
            // rdft overwites the input array ooura_fft_buf with resulting compex number pair for every frequency
            // up to Nyquist frequency which is half the length of the audio. so the output is same size as input.
            rdft(WINDOW, 1, ooura_fft_buf, ooura_ip, ooura_w); 



            double *p = ooura_fft_buf-1;

            for (int j=0; j< image_width; j++) {
                double x = hypot(*++p, *++p);
                // usually spectrograms are log-power-spectrum: log(r**2 + i**2)
                // however, hypot() is documented to be faster than r**2 + i**2
                // and since log(X**2)=2*log(X) and log(X**0.5)=0.5*log(X)
                // and later we normalize all the resulting log() values
                // 2x or 0.5x factor flatens out. Therefore better use the faster hypot()
                x = (x == 0)? -200: log(x); // avoid log(0) division by zero

                // DC frequency up to~25Hz usually don't have meaningful visual audiable information
                // but they may have high energy and may cause reduction of details for higher and more interesting frequencies
                // Here we only ignore DC values -- this may be a small value hence not needed because it is the sum of the normalized signal [-1..+1]
                if (j==0) {
                    x=0; // this is probably not needed since audio signal is normalized, DC should be 0
                    //continue;
                }
                if (x > maxfft) {
                    maxfft = x;
                    //printf("fft power max=%f\n", maxfft);
                } else if (x < minfft) {
                    minfft = x;
                    //printf("fft power min=%f\n", minfft);
                }

                *img++ = (float)x;

            }
        }
        float delfft=maxfft-minfft; // delta
        #ifdef DEBUG
        fprintf(stderr, "TOTAL LINES %d maxfft=%f minfft=%f delfft=%f\n", image_lines, maxfft, minfft, delfft);
        #endif
        img = image_data;

        // now that we have min and max, normalize the spectrogram
        for (int j=0, k=image_width*image_lines; j<k;j++) {
            // *img = *img>maxfft?maxfft:*img // if maxfft ignores DC values, then DC line should be set to zero. otherwise do this histeresis

            float x = *img;
            x = x>maxfft?maxfft:x; // hysterisis is possible in case we skipped some frequences when calculating min and max
            x = (x-minfft)/delfft; // normalized to 0..1 on log scale
            if (x > 1) { // sanity check
                #ifdef DEBUG
                fprintf(stderr, "wrong value %f\n",x);
                #endif
            }
            *img++ = x;
        }
    }

    void spectrogram_histogram_equalization() {
        int pixels = image_width*image_lines;
        int hist[65536] = {0};
        float hist_new[65536] = {0};
        float *img = image_data; // values are normalized float 0..1
        // calculating histogram with 1<<16-1 buckets
        for (int j=0; j<pixels;j++) {
            hist[(int)(65535*(*img++))]++;
        }
        // calculating cumulative frequency and new gray levels
        for (int j = 0, curr = 0; j < 65536; j++) {
            // cumulative frequency
            curr += hist[j];

            // calculating new gray level after multiplying by
            // maximum gray count which is 65535 and dividing by
            // total number of pixels
            hist_new[j] = float(curr) / pixels;
        }
        // performing histogram equalisation by mapping new gray levels
        img = image_data;
        for (int j=0; j<pixels;j++) {
            *img = hist_new[int(65535*(*img))];
            img++;
        }
    }

    void mel_spectrogram(int MELSIZE, int SAMPLE_RATE, float MIN_FREQ, float MAX_FREQ, float MEL_700, float MEL_2595) {

        mel_image_data = new float[image_lines*MELSIZE];
        mel_image_width = MELSIZE;
        auto mel = Mel(image_width, SAMPLE_RATE, MELSIZE, MIN_FREQ, MAX_FREQ, MEL_700, MEL_2595);
        for (int i=0; i<image_lines; i++) {
            mel.transform(image_data+i*image_width, mel_image_data+i*MELSIZE);
        }
    }


    void color(COLORMAPS_t colormap, bool fft_or_mel) {
        Colormaps_t *cm = (Colormaps_t*)COLORMAP[colormap];
        float *img;
        int width;
        if (fft_or_mel) {
            img = image_data;
            width = image_width;
        } else {
            img = mel_image_data;
            width = mel_image_width;
        }
        if (rgb_image_data) delete rgb_image_data;
        rgb_image_data = new float[image_lines*width*3];
        float *rgb = rgb_image_data;

        for (int i=0; i<image_lines; i++) {
            for (int j=0;j<width;j++) {
                int bin_edges = COLORMAP_SIZE-1;
                float x = *img++;
                float seg = 1.0/bin_edges; // 256 equal segments
                float mod_right = (float)fmod(x, seg)/seg; // fraction of a segment
                float mod_left = 1-mod_right;
                int bin_left = bin_edges*x;
                int bin_right = bin_left+1;
                if (bin_left < 0) {
                    #ifdef DEBUG
                    fprintf (stderr, "problem bin_left=%d\n", bin_left);
                    #endif
                   bin_left = 0;
                    bin_right = 1;
                    mod_left=0;
                    mod_right=1;
                 }
                double *l = (*cm)[bin_left];
                double *r = (*cm)[bin_right];
                *rgb++ = *l++*mod_left + *r++*mod_right;
                *rgb++ = *l++*mod_left + *r++*mod_right;
                *rgb++ = *l++*mod_left + *r++*mod_right;
                //*rgb++ = (*cm)[bin_left][0]*mod_left+(*cm)[bin_right][0]*mod_right;
                //*rgb++ = (*cm)[bin_left][1]*mod_left+(*cm)[bin_right][1]*mod_right;
                //*rgb++ = (*cm)[bin_left][2]*mod_left+(*cm)[bin_right][2]*mod_right;
            }
        }
    }
};


