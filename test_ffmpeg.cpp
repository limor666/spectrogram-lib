#include <cstdio>
#include "spectrogram.cpp"

// g++ -DTEST -O6  test_ffmpeg.cpp ../OouraFFT/fftsg.c -lm
// ./a.out guitar.ogg asdf.png
// creates 3 files : fft_asdf.png  fftg_asdf.png  mel_asdf.png
// Requires ffmpeg to be installed. To view the big files use gimp
// Launchs FFmpeg several times, one to read the original audio
// file and another 2 to write the images. In each case,
// data passes between this program and FFmpeg through a pipe.
// 2 colored images are output one for fft and one for mel

#include <chrono> 
int main(int ac, char **av)
{

    if (ac<3) {
	    printf("usage: %s audiofile.ogg filename.png\n", av[0]);
	    exit(1);
    }
    // expect input in ogg format
    // output 1 channel float 32 bit int little-endian
    FILE *pipe;
    char str[200];
    // test: ffmpeg -loglevel warning -acodec vorbis -i guitar.ogg -f f32le -ar 44100 -ac 1 - | play -t raw -r 44100 -e float -b 32 -c 1 -
    sprintf(str, "ffmpeg -loglevel warning -acodec vorbis -i %s -f f32le -ar 44100 -ac 1 -", av[1]);
    printf("Input command: %s\n", str);
    pipe  = popen(str, "r");

    // read the piped audio file into an audio buffer
    // resize the buffer every 1000000 bytes read
    int bufsize=1000000;
    int buffers=1;
    float *audio = (float*)malloc(bufsize);
    int index = 0;
    while (!feof(pipe)) {
        int c = fread((char*)audio+index, 1, bufsize-(index%bufsize), pipe);
        if (!c) continue;
        index += c;
        if (index%bufsize == 0) {
            audio = (float*)realloc(audio, ++buffers*bufsize);
        }
    }
    pclose(pipe);
    fprintf(stderr, "Audio samples: %d, leftover bytes: %d\n", index/(int)sizeof(float), index%(int)sizeof(float));
    int samples = index/4;

    #define NOW() std::chrono::high_resolution_clock::now() // c++ sucks
    #define DURATION(end,start) std::chrono::duration<double, std::milli>(end-start).count()/1000 // c++ sucks
    #define PRINTDURATION(str) time2=NOW(); fprintf(stderr, "%-22s %4.2f sec\n", str, DURATION(time2, time1)); time1 = time2;
    auto time1 = NOW();
    auto time2 = NOW();

    // Create the Spectrogram class instance
    auto spec = spectrogram(audio, samples);
    PRINTDURATION("spectrogram class init");
    spec.normalize_audio();
    PRINTDURATION("normalize");
    spec.spectrogram_from_audio(-512, 250, 1024, fft_window::HAMMING);
    PRINTDURATION("spectrogram");
    spec.spectrogram_histogram_equalization();
    PRINTDURATION("histogram equalization");

    // Output GRAY fft image via a pipe to ffmpeg which accepts incantations to generate an image
    sprintf(str, "ffmpeg -loglevel warning -y -f rawvideo -s %dx%d -pix_fmt gray16le -i - fftg_%s", spec.image_width, spec.image_lines, av[2]);
    printf("Output command: %s\n", str);
    pipe = popen(str, "w");
    for (struct { int i; float *u; uint16_t x; } p = { 0, spec.image_data, 0 }; p.i<spec.image_lines*spec.image_width; p.i++) {
        p.x = (uint16_t)(65535*(*p.u++));
        fwrite(&p.x, 1, sizeof(uint16_t), pipe); // not sure how efficient it is to call fwrite a million times
    }
    pclose(pipe);

    PRINTDURATION("wrote gray image");
    // add some color. create the color image buffer
    spec.color(COLORMAP_MAGMA, true);
    PRINTDURATION("color");

    // Output RGB fft image via a pipe to ffmpeg which accepts incantations to generate an image
    sprintf(str, "ffmpeg -loglevel warning -y -f rawvideo -s %dx%d -pix_fmt rgb48le -i - fft_%s", spec.image_width, spec.image_lines, av[2]);
    printf("Output command: %s\n", str);
    pipe = popen(str, "w");
    for (struct { int i; float *u; uint16_t x; } p = { 0, spec.rgb_image_data, 0 }; p.i<spec.image_lines*spec.image_width*3; p.i++) {
        p.x = (uint16_t)(65535*(*p.u++));
        fwrite(&p.x, 1, sizeof(uint16_t), pipe); // not sure how efficient it is to call fwrite a million times
    }
    pclose(pipe);

    PRINTDURATION("wrote color image");
    // distort the GRAY image with mel algorithm
    spec.mel_spectrogram(500, 44100, 30.0, 0.7*44100/2, 1000.0, 3000.0);
    PRINTDURATION("mel");
    spec.color(COLORMAP_MAGMA, false);
    PRINTDURATION("color");


    // Output mel image via a pipe to ffmpeg which accepts incantations to generate an image
    sprintf(str, "ffmpeg -loglevel warning -y -f rawvideo -s %dx%d -pix_fmt rgb48le -i - mel_%s", spec.mel_image_width, spec.image_lines, av[2]);
    printf("Output command: %s\n", str);
    pipe = popen(str, "w");
    for (struct { int i; float *u; uint16_t x; } p = { 0, spec.rgb_image_data, 0 }; p.i<spec.image_lines*spec.mel_image_width*3; p.i++) {
        p.x = (uint16_t)(65535*(*p.u++));
        fwrite(&p.x, 1, sizeof(uint16_t), pipe); // not sure how efficient it is to call fwrite a million times
    }
    pclose(pipe);
    PRINTDURATION("wrote color image");


}
