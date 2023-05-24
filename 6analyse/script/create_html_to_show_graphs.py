from pathlib import Path
from string import Template


class HtmlTemplate(Template):
    delimiter = "~"


def create_html(paths, output):
    img_template = Template(
        '<li><img src="$url?width=300&height=300" alt="$alt" /></li>'
    )
    with open("script/template.html") as f:
        html_template = HtmlTemplate(f.read())
    img_tags = []
    for p in paths:
        img_tags.append(
            img_template.substitute(alt=p, url="/resized" / Path(*p.parts[1:]))
        )

    with open(output, mode="w", encoding="utf8") as f:
        f.write(html_template.substitute(images="".join(img_tags)))


files = [Path(*Path(f).parts[1:]) for f in snakemake.input]  # type: ignore
output = snakemake.output[0]  # type: ignore
create_html(files, output)
