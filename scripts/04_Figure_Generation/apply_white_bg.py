from PIL import Image

out_dir = r"E:\Zika\results\phaseF\Drug_Repurposing_Results"
png_file = rf"{out_dir}\Figure15_Sankey_HighRes.png"

# Open the image
img = Image.open(png_file).convert("RGBA")

# Create a white background image
bg = Image.new("RGBA", img.size, "WHITE")
bg.paste(img, (0, 0), img)

# Convert to RGB (dropping alpha)
rgb_img = bg.convert("RGB")

# Save as PNG with 600 DPI
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.png", dpi=(600, 600))

# Save as JPG and TIFF with 600 DPI
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.jpg", "JPEG", quality=100, dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.tif", "TIFF", dpi=(600, 600), compression="tiff_lzw")
print("Done!")
