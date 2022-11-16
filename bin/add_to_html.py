#!/usr/bin/env python3
# encoding: utf-8
import base64


def get_base64_encoded_image(image_path):
    with open(image_path, "rb") as img_file:
        return 'data:image/png;base64,' + base64.b64encode(img_file.read()).decode('utf-8')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--image", help="Path to image to append to html", action='append')
    parser.add_argument("--html", help="Path to HTML file")
    parser.add_argument("--out", help="Path to output HTML file")

    args = parser.parse_args()
    images = args.image
    html = args.html
    out = args.out
    with open(html, 'rt') as f:
        html_text = f.read()

    text_to_add = '<div class="container-fluid">'
    for i in range(len(images)):
        image = images[i]
        # title = titles[i]
        # if title != 'na':
        #     text_to_add += '<h4>{}</h4>'.format(title)
        b64 = get_base64_encoded_image(image)
        text_to_add += '<img src="{}">'.format(b64)
        text_to_add += '<br />'
    text_to_add += '</div>'
    index = html_text.rfind('</body>')
    html_text = html_text[:index] + text_to_add + html_text[index:]
    with open(out, 'wt') as f:
        f.write(html_text)
