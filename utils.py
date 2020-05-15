from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
from PIL import Image
import cairosvg
import base64
import os
from io import BytesIO
from tqdm import tqdm_notebook

def rgb2rgba(image):
    img = image.convert("RGBA")
    datas = img.getdata()
    newData = []
    for item in datas:
        if item[0] == 255 and item[1] == 255 and item[2] == 255:
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)
    img.putdata(newData)
    return img
    
def load_pics(smiles_list):
    DrawingOptions.dotsPerAngstrom = 70
    DrawingOptions.bondLineWidth = 2
    DrawingOptions.atomLabelFontSize = 28

    pic_urls = []
    print('Loading molecular structures...')
    for smiles in tqdm_notebook(smiles_list):
        try:
            Draw.MolToFile(Chem.MolFromSmiles(smiles), "temp.png", size=(160, 100))
            image = rgb2rgba(Image.open('temp.png'))

            buffered = BytesIO()
            image.save(buffered, format="png")       
            url = 'data:image/png;base64,' + base64.b64encode(buffered.getvalue()).decode('utf-8')
            pic_urls.append(url)

            os.remove("temp.svg")
            os.remove("temp.png")
        except:
            pic_urls.append(smiles)
    return pic_urls

