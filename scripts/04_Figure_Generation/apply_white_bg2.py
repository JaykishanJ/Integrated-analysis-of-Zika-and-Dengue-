from PIL import Image

out_dir = r"E:\Zika\results\phaseF\Drug_Repurposing_Results"
png_file = rf"{out_dir}\Figure15_Sankey.png"

# Open the image
img = Image.open(png_file).convert("RGBA")

# Create a white background image
bg = Image.new("RGBA", img.size, "WHITE")
bg.paste(img, (0, 0), img)

# Convert to RGB (dropping alpha)
rgb_img = bg.convert("RGB")

# Save as PNG with 600 DPI (overwriting original)
rgb_img.save(png_file, dpi=(600, 600))
print("Overwrote original PNG")
