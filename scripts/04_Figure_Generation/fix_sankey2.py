import re

filepath = r'e:\Zika\results\phaseF\Drug_Repurposing_Results\PhaseF_Drug_repurpusing.R'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Remove Thioridazine lines completely
content = re.sub(r'\s*"Thioridazine" = "Thioridazine",?', '', content, flags=re.IGNORECASE)
content = re.sub(r'\s*"Thioridazine",?', '', content, flags=re.IGNORECASE)
content = re.sub(r'data\.frame\([^)]*Thioridazine[^)]*\)', lambda m: m.group(0).replace(', "Thioridazine"', '').replace('"Thioridazine",', ''), content)

# Remove the specific edge for Thioridazine in sankey data if it exists
content = re.sub(r'c\([^)]*Thioridazine[^)]*\)', lambda m: m.group(0).replace(', "Thioridazine"', '').replace('"Thioridazine",', ''), content)

# Change NAC to N-acetylcysteine
content = content.replace('"N-Acetylcysteine" = "NAC"', '"N-Acetylcysteine" = "N-acetylcysteine"')

# Fix webshot zoom to 6 and add white background / 600 DPI as before
old_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600)
      cat("Static images (PNG and PDF) saved using", pkg, "\\n")'''
      
new_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600, zoom = 6)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600, zoom = 6)
      
      # Python script will be called externally to process the image to avoid magick issues.
'''
content = content.replace(old_webshot, new_webshot)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(content)

print("Modified script.")
