#include <cmath>

class fft_window {

public:
    typedef enum {
        HANNING,
        HAMMING,
        TRIANGULAR,
        GAUSS,
        BLACKMAN_HARRIS,
        FLAT,
        RANDOM,
        ALMOST_FLAT
    } Type; 

    double *data;

    fft_window(int length, Type type) {

        data = new double[length](); // set to zero

        switch (type) {

            case HANNING:
                for (int i=0; i< length; i++)
                    data[i] = Window_hanningAtPoint(i, length);
                break;

            case HAMMING:
                for (int i=0; i< length; i++)
                    data[i] = Window_hammingAtPoint(i, length);
                break;

            case TRIANGULAR:
                for (int i=0; i< length; i++)
                    data[i] = Window_triangularAtPoint(i, length);
                break;

            case GAUSS:
                for (int i=0; i< length; i++)
                    data[i] = Window_gaussAtPoint(i, length);
                break;

            case BLACKMAN_HARRIS:
                for (int i=0; i< length; i++)
                    data[i] = Window_blackmanHarrisAtPoint(i, length);
                break;

            case FLAT:
                for (int i=0; i< length; i++)
                    data[i] = 1.0f;
                break;

            case RANDOM:
                for (int i=0; i< length; i++)
                    data[i] = Window_randomAtPoint(i, length);
                break;

            case ALMOST_FLAT:
                for (int i=0; i< length; i++){
                    data[i] =  powf( 2 * Window_hanningAtPoint(i, length), 0.3 );
                    data[i] = powf( data[i] , 1.5); 
                    if ( data[i] > 1.0f ) data[i] = 1.0f;
                }
                break;
                ;
            default:
                for (int i=0; i< length; i++)
                    data[i] = 1.0f;
                //printf("Undefined dataType!\n");
                break;
        }
    }

private:


    float Window_hanningAtPoint(float point, int numSamples){
        return 0.5f * ( 1.0f - cosf ( (2.0f * M_PI * point) / ( (float)numSamples - 1.0f) ) );
    }


    float Window_hammingAtPoint(int point, int numSamples){
        return 0.54f - 0.46f * cosf ( (2.0f * M_PI * point) / ( (float)numSamples - 1.0f) ) ;
    }


    float Window_triangularAtPoint(int point, int numSamples){
        return ( 2.0f / numSamples ) * ( ( numSamples * 0.5f ) - fabs( point - ( numSamples -1.0f ) * 0.5f ) );
    }


    float Window_gaussAtPoint(int point, int numSamples){
        
        float bellW = 0.4f;
        return powf ( M_E, -0.5f * powf ( ( point - ( numSamples - 1 ) * 0.5f ) / ( bellW * ( numSamples - 1 ) * 0.5f ) , 2.0f ) );
    }


    float Window_blackmanHarrisAtPoint(int point, int numSamples){
        
        return 0.35875f		- 0.48829f * cosf( 2.0f * M_PI * point / (numSamples-1) ) 
                            + 0.14128f * cosf( 4.0f * M_PI * point / (numSamples-1) ) 
                            - 0.01168f * cosf( 6.0f * M_PI * point / (numSamples-1) );
    }


    float Window_randomAtPoint(int point, int numSamples){
        return rand()%1000 / 1000.0f;
    }
};
