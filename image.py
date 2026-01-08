#!/usr/bin/python


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import oskar

base = '/gpfs01/home/ppxjf3/OSKAR/'

image_type = sys.argv[1]
freq = float(sys.argv[2]) #MHz 
sky_root_name = sys.argv[3]
fov_deg  = float(sys.argv[4])#5.0
num_pixels_side = int(sys.argv[5])#512
field = sys.argv[7]
telescope = sys.argv[6]
telescope_model = sys.argv[8]
weighting = sys.argv[9]

imager = oskar.Imager()
ms=".vis"
ms_name = base + 'vis/%s_%s_%s_%s_%03.1f_%03d.vis' % (sky_root_name, telescope, telescope_model, field, fov_deg, num_pixels_side)
out_name = base + '%s_%s_%s_%s_%03.1f_%04d_%s_%s'% (sky_root_name, telescope, telescope_model, field, fov_deg, num_pixels_side,weighting,image_type)
imager.set(image_type=image_type, fov_deg=fov_deg, image_size=num_pixels_side)
imager.set(input_file=ms_name, output_root=out_name)
output = imager.run(return_images=1)
image = output["images"][0]

# Render the image using matplotlib and save it as a PNG file.
im = plt.imshow(image, cmap="jet")
plt.gca().invert_yaxis()
plt.colorbar(im)
plt.savefig(out_name+'.png')
