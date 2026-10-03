import sys

filepath = r'e:\Zika\results\phaseF\Drug_Repurposing_Results\PhaseF_Drug_repurpusing.R'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Change NAC to N-acetylcysteine
content = content.replace('"N-Acetylcysteine" = "NAC"', '"N-Acetylcysteine" = "N-acetylcysteine"')

# Replace thioridazine in drug_display list
content = content.replace('"Disulfiram" = "Disulfiram",\n  "Thioridazine" = "Thioridazine",\n  "N-Acetylcysteine" = "N-acetylcysteine"', '"Disulfiram" = "Disulfiram",\n  "N-Acetylcysteine" = "N-acetylcysteine"')
content = content.replace('"Thioridazine" = "Thioridazine",\n  "N-Acetylcysteine" = "N-acetylcysteine"', '"N-Acetylcysteine" = "N-acetylcysteine"')
content = content.replace('"Thioridazine" = "Thioridazine",', '')
content = content.replace('"Thioridazine",', '')

# Replace in drug_score_data
content = content.replace('"Tamoxifen", "Disulfiram", "Thioridazine"),', '"Tamoxifen", "Disulfiram"),')
content = content.replace('Evidence_Score = c(6, 5, 5, 5, 5, 5, 3, 5, 3, 1),', 'Evidence_Score = c(6, 5, 5, 5, 5, 5, 3, 5, 3),')

# Remove from drug vector if exists
content = content.replace('c("Methotrexate", "Bortezomib", "Lithium Carbonate", "Metformin",\n                        "Acadesine", "N-Acetylcysteine", "Midostaurin",\n                        "Tamoxifen", "Disulfiram", "Thioridazine")', 'c("Methotrexate", "Bortezomib", "Lithium Carbonate", "Metformin",\n                        "Acadesine", "N-Acetylcysteine", "Midostaurin",\n                        "Tamoxifen", "Disulfiram")')

# Replace webshot
old_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600)'''
new_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600, zoom=6)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600, zoom=6)'''
content = content.replace(old_webshot, new_webshot)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(content)
print("done editing")
