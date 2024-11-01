// Thank you chatgpt

#include <cmath>
#include <cstdio>
//#include <stdio.h>
//#define _GNU_SOURCE
//#include <math.h>
//#include <stdlib.h>

// https://en.wikipedia.org/wiki/Mel_scale
// The mel concept: log10(1.0->2.4) = 0->0.4 is somewhat linear and starts acting "log" thereafter
// therefore log10(1.0+freq/700) is somewhat linear up to freq=1000 
// 2595*log10(1+freq/700) ~= 1000 when freq=1000
// mel value for freq 0 -> 1000 are roughly like freq
// mel value for freq 1000+ become log squashed
// for freq=22000 (nyquist), mel = 3920
// mel transform of FFT vector, magnifies freq 0..1000 over 50% of the resulting vector values and freq 1000..22000 are log squashed over remaining 50%
// For the purpose of voice visualization this is probably good. for some drums better use full spectrum FFT.
// However, in the API we add the possibility of setting the min/max frequencies and linear cutoff
//
// The bigger the MEL_700 the less smeared the lower frequencies
// The value of MEL_2595 doesn't have any effect at all. don't have patience to explore why




class Mel 
{

private:
    float *filterbank;
    int *bin_points; // this is global to be used in apply_filterbank()
    int FFT_SIZE;  // Size of the FFT (e.g., 512, 1024)
    int SAMPLE_RATE;  // Sampling rate of the signal (e.g., 16000Hz)
    int NUM_MEL_BINS; // Number of Mel bins (mel output )
    float MIN_FREQ=50; // ignore FFT values less than this frequency
    float MAX_FREQ=10000; // ignore FFT values more than this frequency
    float MEL_700=700;
    float MEL_2595=2595;  // seems to have no effect




// Function to convert a frequency to the Mel scale
float freq_to_mel(float freq) {
    return MEL_2595 * log10(1.0 + freq / MEL_700);
    //return 1127.0 * log10(1.0 + freq / 700.0);
    //return 5000 * log10(1.0 + freq / 1500.0);
}

// Function to convert from Mel scale to frequency
float mel_to_freq(float mel) {
    return MEL_700 * (exp10(mel / MEL_2595) - 1.0);
    //return 700.0 * (exp10(mel / 1127) - 1.0);
    //return 1500 * (exp10(mel / 5000) - 1.0);
}

// Function to compute the Mel filterbank
// It maps FFT bins to the mel scale and sets up triangular filters to sum
// the energies in the FFT bins that correspond to each mel bin.
void create_filterbank() {
    int i, j;

    // Frequency range of FFT
    float MIN_MEL = freq_to_mel(MIN_FREQ);  // Minimum mel
    float max_mel = freq_to_mel(MAX_FREQ);  // Nyquist frequency in mel

    // Create mel points spaced linearly between MIN_FREQ and MAX_FREQ
    float mel_points[NUM_MEL_BINS + 2];
    for (i = 0; i < NUM_MEL_BINS + 2; i++) {
        mel_points[i] = MIN_MEL + i * (max_mel - MIN_MEL) / (NUM_MEL_BINS + 1);
    }

    // Convert mel points to frequency points
    float freq_points[NUM_MEL_BINS + 2];
    for (i = 0; i < NUM_MEL_BINS + 2; i++) {
        freq_points[i] = mel_to_freq(mel_points[i]);
    }

    // Convert frequency points to FFT bin numbers
    for (i = 0; i < NUM_MEL_BINS + 2; i++) {
        bin_points[i] = (int)(FFT_SIZE * freq_points[i] / (SAMPLE_RATE/2));
    }

    // Create the filterbank
    for (i = 0; i < NUM_MEL_BINS; i++) {
        #ifdef DEBUG
        printf("filter %3d: %3d %3d -- %7.1f %7.1f [", i, bin_points[i], bin_points[i+2], freq_points[i], freq_points[i+2]);
        #endif
        int asdf = bin_points[i+2] - bin_points[i];
        switch(asdf) {
            case 0:
                filterbank[i * FFT_SIZE + bin_points[i]] = 1;
                #ifdef DEBUG
                printf("%3d:1.0 ", bin_points[i]);
                #endif
                break;
            /*case 1: 
                filterbank[i * FFT_SIZE + bin_points[i]+1] = 1;
                printf("%3d:1.0 ", bin_points[i]+1);
                break;
                */
            default:
                float mm = 0;
                for (j = bin_points[i]; j <= bin_points[i+2]; j++) {
                    //float m = sin(M_PI*(j-bin_points[i]+1)/(asdf+2));
                    float m = 0.5-fabs(j-bin_points[i]+1-(asdf+2)/2.0)/(asdf+2);
                    mm += m;
                    filterbank[i*FFT_SIZE+j] = m;
                }
                for (j = bin_points[i]; j <= bin_points[i+2]; j++) { // normalize the values so that the sum is 1
                    filterbank[i*FFT_SIZE+j] /= mm;
                    #ifdef DEBUG
                    printf("%3d:%3.2f ", j, filterbank[i * FFT_SIZE + j]);
                    #endif
                }
	    }
        #ifdef DEBUG
	    printf("]\n");
        #endif
    }
}




public:
    Mel(
        int FFT_SIZE_,  // Size of the FFT (e.g., 512, 1024)
        int SAMPLE_RATE_,  // Sampling rate of the signal (e.g., 16000Hz)
        int NUM_MEL_BINS_, // Number of Mel bins (mel output )
        float MIN_FREQ_=50, // ignore FFT values less than this frequency
        float MAX_FREQ_=10000, // ignore FFT values more than this frequency
        float MEL_700_=700,
        float MEL_2595_=2595)  // seems to have no effect
    {
        FFT_SIZE = FFT_SIZE_;
        SAMPLE_RATE = SAMPLE_RATE_;
        NUM_MEL_BINS = NUM_MEL_BINS_;
        MIN_FREQ = MIN_FREQ_;
        MAX_FREQ = MAX_FREQ_;
        MEL_700 = MEL_700_;
        MEL_2595 = MEL_2595_;
        filterbank = new float[NUM_MEL_BINS * FFT_SIZE];
        bin_points = new int[NUM_MEL_BINS * NUM_MEL_BINS+2];
        create_filterbank();
    }
    ~Mel() {
        delete[] filterbank;
        delete[] bin_points;
    }


    // Function to apply Mel filterbank to FFT vector
    // Generates the mel-frequency vector.
    // After applying the filterbank, the function computes the log energy of each mel bin.
    // before this optimization, it did matrix multiplication: mel_output[i] += fft_output[j] * filterbank[i * FFT_SIZE + j]; 
    // after optimization we get 2x improvements. Must use: gcc -O6
    void transform(float* fft_power, float* mel) {
        int i, j;

        for (i = 0; i < NUM_MEL_BINS; i++) {
            mel[i] = 0.0;
            int i1 = bin_points[i]; // first one is always 0
            int i2 = bin_points[i+2];
            float *fb = &filterbank[i*FFT_SIZE+i1];
            float *ff = &fft_power[i1];
            for (j = i1; j <= i2; j++) {
                mel[i] += (*(ff++)) * (*(fb++));
            }
        }
    }

};

