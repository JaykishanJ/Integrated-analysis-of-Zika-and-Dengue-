from PIL import Image, ImageDraw, ImageFont

out_dir = r"E:\Zika\results\phaseF\Drug_Repurposing_Results"
png_file = rf"{out_dir}\Figure15_Sankey.png"

# Open the image
img = Image.open(png_file).convert("RGBA")

# We want to add a title. Create a new taller image.
# Image is 1800 x 1600. Title area = 150px height.
# And we also want to crop out the top 50px if there's any stray pixels, but since we removed the text, it should be clean.
new_height = img.height + 150
bg = Image.new("RGBA", (img.width, new_height), "WHITE")

# Paste original image lower
bg.paste(img, (0, 150), img)

# Draw title
draw = ImageDraw.Draw(bg)
# Using default font because we don't know Arial location easily, but we can make it large by using ImageFont.truetype if available
try:
    font = ImageFont.truetype("arial.ttf", 60)
except:
    font = ImageFont.load_default()

title_text = "Hub Gene-Drug Candidate Network"
# calculate text size
bbox = draw.textbbox((0, 0), title_text, font=font)
text_width = bbox[2] - bbox[0]
text_height = bbox[3] - bbox[1]

x = (img.width - text_width) / 2
y = 60

draw.text((x, y), title_text, fill="black", font=font)

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

print("Saved all formats successfully with drawn title.")
