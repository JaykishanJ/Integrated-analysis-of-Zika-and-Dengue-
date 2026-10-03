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

# Save all the original filenames with 600 DPI
rgb_img.save(rf"{out_dir}\Figure15_Sankey.png", dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey.jpg", "JPEG", quality=100, dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey.tif", "TIFF", dpi=(600, 600), compression="tiff_lzw")

# Also save the _HighRes filenames just in case
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.png", dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.jpg", "JPEG", quality=100, dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.tif", "TIFF", dpi=(600, 600), compression="tiff_lzw")

print("Saved all formats successfully with title.")
