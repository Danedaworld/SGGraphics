
![snowman](https://raw.githubusercontent.com/Danedaworld/SGGraphics/master/Other%20Images/reflective%20snowman.bmp)
![reflective walls](https://raw.githubusercontent.com/Danedaworld/SGGraphics/master/Portfolio/Indirect_Reflective_Wall%201000%20Samples.bmp)

5.1 - OpenMP support enabled. I use each thread to process a line of pixels.

5.2 - Antialiasing super-samples each pixel with 9 equidistant, precomputed offsets. I included a simple image that demonstrates this. It uses the vanilla config file provided in the sample.

5.3 - Area light and Soft Shadows - IMPORTANT: I changed EMIT to actually take in arbitrary float values. The reason for this is so I can use it as a emmitance value. I send about 30 randomly
generated shadow-feelers and then take a illuminance weight to get the expected soft shadow result. I have included a simple image that demonstrates this. It uses the vanilla config file provided.

6 - MONTE CARLO PATH TRACING - WITH DIRECT RAY TRACING -

I use Monte-Carlo path tracing to render some images. The output image is stored as Indirect.bmp Direct.bmp and Combined.bmp. Indirect is as you'd expect, vanilla MC path tracing. On Line 562, if you change 'hitcount < -1' to 'hitcount < 2' it will ONLY obtain indirect color based on what was stated on the piazza posts. The Direct.bmp is a vanilla ray trace. I use a large cache to store all the light-sampled weights for combining Indirect and Direct lights. Since light-sampling will yield a LOT of 100% sampling for perfectly illuminated place, I bias it towards 50%. Using this biasing, I can get some indirect lighting to show up in the direct + indirect combination. The output of this file is Combination.bmp. I have included some example images.

Note that you can set floating point values for MIRR and TRANS to get weightings for them but this will only work with vanilla MC path tracing. I have included a cool example.


Timings:
Without OpenMP - 10.35 seconds
With 8-cores - 2.89 seconds

Averaged over 5 trials