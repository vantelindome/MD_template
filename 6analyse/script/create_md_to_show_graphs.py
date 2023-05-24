from pathlib import Path


def create_markdown(files, markdown_file):
    img_tags = [f"## [{file}]({file})\n![]({file})" for file in sorted(files)]
    markdown = "\n".join(img_tags)
    with open(markdown_file, "w") as f:
        f.write(markdown)


files = [Path(*Path(f).parts[1:]) for f in snakemake.input]  # type: ignore
output = snakemake.output[0]  # type: ignore
create_markdown(files, output)
