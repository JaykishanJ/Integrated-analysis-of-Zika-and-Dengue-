import re

def fix_r_script(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except:
        return
        
    old_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600)
      cat("Static images (PNG and PDF) saved using", pkg, "\\n")'''
      
    new_webshot = '''getExportedValue(pkg, "webshot")(html_file, file = png_file, vwidth = 1800, vheight = 1600, zoom = 6)
      getExportedValue(pkg, "webshot")(html_file, file = pdf_file, vwidth = 1800, vheight = 1600, zoom = 6)
      
      # Convert to JPG and TIF using magick if available
      if (requireNamespace("magick", quietly = TRUE)) {
          img <- magick::image_read(png_file)
          magick::image_write(img, path = file.path(out_dir, "Figure15_Sankey.jpg"), format = "jpeg", quality = 100)
          magick::image_write(img, path = file.path(out_dir, "Figure15_Sankey.tiff"), format = "tiff", compression = "LZW")
          cat("Saved high-res PNG, PDF, JPG, and TIFF (approx 600 DPI).\\n")
      } else {
          cat("Saved high-res PNG and PDF. Install 'magick' for JPG and TIFF formats.\\n")
      }
'''
    
    content = content.replace(old_webshot, new_webshot)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

fix_r_script(r'E:\Zika\scripts\PhaseF_1_Drug_repurpusing.R')
print("done")
