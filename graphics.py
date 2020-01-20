from Lib.PIL import Image
import solve
import gauss

def save_image(image, path):
  image.save(path, 'png')

def create_image(i, j):
  image = Image.new("RGB", (i, j), "white")
  return image

def to_gray_pixel(value, max_value):
  red_value = int(255 * value / max_value)
  return (255 - red_value, 255 - red_value, 255 - red_value)

def graphical_solve(probes, img_name):
  image = create_image(probes*2, probes*2)
  pixels = image.load()
  values = solve.ex_u_values(probes)
  max_value = max([max(row, key=lambda v: v if v is not None else float('-inf')) for row in values])
  for (i, row) in enumerate(values):
    for (j, value) in enumerate(row):
      if value is not None: pixels[j, i] = to_gray_pixel(value, max_value)
  save_image(image, img_name)

graphical_solve(320, 'rozw.png')
