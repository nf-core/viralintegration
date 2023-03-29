#!/usr/bin/env python3
# encoding: utf-8

#########################################
# Author: [brianjohnhaas](https://github.com/brianjohnhaas)
# File: add_to_html.py
# Source: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/master/util/add_to_html.py
# Source+commit: https://github.com/broadinstitute/CTAT-VirusIntegrationFinder/blob/3082b0a027e5f6c224cc9fd3cdbe8bf72fb04ec7/util/add_to_html.py
# Download Date: 2022-12-28, commit: 3082b0a
# This source code is licensed under the BSD 3-Clause license
#########################################

import base64


def get_base64_encoded_image(image_path):
    with open(image_path, "rb") as img_file:
        return "data:image/png;base64," + base64.b64encode(img_file.read()).decode("utf-8")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--image", help="Path to image to append to html", action="append")
    parser.add_argument("--html", help="Path to HTML file")
    parser.add_argument("--out", help="Path to output HTML file")

    args = parser.parse_args()
    images = args.image
    html = args.html
    out = args.out
    with open(html, "rt") as f:
        html_text = f.read()

    text_to_add = '<div class="container-fluid">'
    for i in range(len(images)):
        image = images[i]
        # title = titles[i]
        # if title != 'na':
        #     text_to_add += '<h4>{}</h4>'.format(title)
        b64 = get_base64_encoded_image(image)
        text_to_add += '<img src="{}">'.format(b64)
        text_to_add += "<br />"
    text_to_add += "</div>"
    index = html_text.rfind("</body>")
    html_text = html_text[:index] + text_to_add + html_text[index:]
    with open(out, "wt") as f:
        f.write(html_text)
