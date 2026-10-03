from PIL import Image

out_dir = r"E:\Zika\results\phaseF\Drug_Repurposing_Results"
files = ["Figure15_Sankey.png", "Figure15_Sankey.jpg", "Figure15_Sankey.tif"]

for f in files:
    img = Image.open(f"{out_dir}\\{f}")
    print(f"{f}: DPI={img.info.get('dpi', 'Not Set')}, Size={img.size}")
