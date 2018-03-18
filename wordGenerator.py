"""
Using a mask, you can generate wordclouds in arbitrary shapes.
(Based off of example files from https://github.com/amueller/word_cloud)
"""

from os import path
from PIL import Image
import numpy as np
import csv

from wordcloud import WordCloud, STOPWORDS

d = path.dirname(__file__)

# Read the whole text.
for i in range(1, 107):
    print(i)
    file = 'tad_GO/tad' + str(i) + '.txt'
    text = ""
    try:
        input_reader = open(path.join(d, file))

        for row in csv.reader(input_reader, delimiter='\t'):
            try:
                text += row[3] + " "
            except IndexError:
                text += row[1] + " "
        input_reader.close()

        # read the mask image
        circle_mask = np.array(Image.open(path.join(d, "Black_Circle.jpg")))

        stopwords = set(STOPWORDS)
        stopwords.add("the")
        stopwords.add("protein")
        stopwords.add("activity")
        stopwords.add("binding")
        stopwords.add("cell")
        stopwords.add("process")

        wc = WordCloud(background_color="white", max_words=2000, mask=circle_mask,
               stopwords=stopwords, prefer_horizontal=1)
        # generate word cloud
        wc.generate(text)
        # store to file
        wc.to_file(path.join(d, "tad-cloud/tad" + str(i) + "-cloud.png"))

    except IOError:
        print(str(i) + " is an empty TAD")
