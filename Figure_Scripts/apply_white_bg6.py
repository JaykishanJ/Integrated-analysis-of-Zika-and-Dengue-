from PIL import Image, ImageDraw, ImageFont

out_dir = r"E:\Zika\results\phaseF\Drug_Repurposing_Results"
png_file = rf"{out_dir}\Figure15_Sankey.png"

img = Image.open(png_file).convert("RGBA")

# Adjust title font size to match the smaller 18px node text scale. 
# Node text 18px at 6x zoom = 108px. Title should be around 140px.
font_size = 140
font_path = r"C:\Windows\Fonts\arialbd.ttf"
try:
    font = ImageFont.truetype(font_path, font_size)
except:
    font = ImageFont.load_default()

top_margin = 300
new_height = img.height + top_margin
bg = Image.new("RGBA", (img.width, new_height), "WHITE")

bg.paste(img, (0, top_margin), img)

draw = ImageDraw.Draw(bg)
title_text = "Hub Gene-Drug Candidate Network"
bbox = draw.textbbox((0, 0), title_text, font=font)
text_width = bbox[2] - bbox[0]
text_height = bbox[3] - bbox[1]

x = (img.width - text_width) / 2
y = (top_margin - text_height) / 2

draw.text((x, y), title_text, fill="black", font=font)

rgb_img = bg.convert("RGB")

rgb_img.save(rf"{out_dir}\Figure15_Sankey.png", dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey.jpg", "JPEG", quality=100, dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey.tif", "TIFF", dpi=(600, 600), compression="tiff_lzw")

rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.png", dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.jpg", "JPEG", quality=100, dpi=(600, 600))
rgb_img.save(rf"{out_dir}\Figure15_Sankey_HighRes.tif", "TIFF", dpi=(600, 600), compression="tiff_lzw")

print("Saved all formats successfully with 18px-equivalent scaled title.")
