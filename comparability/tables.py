import re
import xml.etree.ElementTree as ET


def _strip_ns(elem):

    for e in elem.iter():
        if isinstance(e.tag, str) and "}" in e.tag:
            e.tag = e.tag.split("}", 1)[1]
    return elem


def parse_article(path):

    try:
        tree = ET.parse(path)
    except ET.ParseError:

        with open(path, encoding="utf-8") as fh:
            raw = fh.read()
        raw = re.sub(r"&(?!(?:amp|lt|gt|quot|apos|#\d+|#x[0-9a-fA-F]+);)", "&amp;", raw)
        try:
            tree = ET.ElementTree(ET.fromstring(raw))
        except ET.ParseError:
            return None
    return _strip_ns(tree.getroot())


def cell_text(elem):

    parts = []
    for node in elem.iter():
        if node.tag in ("sub",):
            if node.text:
                parts.append("_" + node.text)
        elif node.tag in ("sup",):
            if node.text:
                parts.append("^" + node.text)
        else:
            if node.text:
                parts.append(node.text)
        if node is not elem and node.tail:
            parts.append(node.tail)
    text = "".join(parts)
    text = text.replace(" ", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _rows_of(table_elem):

    for thead in table_elem.findall(".//thead"):
        for tr in thead.findall(".//tr"):
            yield "head", tr
    for tbody in table_elem.findall(".//tbody"):
        for tr in tbody.findall(".//tr"):
            yield "body", tr

    if table_elem.find(".//thead") is None and table_elem.find(".//tbody") is None:
        for tr in table_elem.findall(".//tr"):
            yield "body", tr


def build_grid(table_elem):

    occupied = {}
    head, body = [], []
    row_idx = 0

    for section, tr in _rows_of(table_elem):
        row = []
        col = 0
        for cell in tr:
            if cell.tag not in ("td", "th"):
                continue

            while (row_idx, col) in occupied:
                row.append(occupied.pop((row_idx, col)))
                col += 1
            text = cell_text(cell)
            try:
                colspan = max(1, int(cell.get("colspan", 1)))
            except ValueError:
                colspan = 1
            try:
                rowspan = max(1, int(cell.get("rowspan", 1)))
            except ValueError:
                rowspan = 1
            for c in range(colspan):
                row.append(text)
                for r in range(1, rowspan):
                    occupied[(row_idx + r, col + c)] = text
            col += colspan

        while (row_idx, col) in occupied:
            row.append(occupied.pop((row_idx, col)))
            col += 1

        (head if section == "head" else body).append(row)
        row_idx += 1

    width = max([len(r) for r in head + body] or [0])
    pad = lambda rows: [r + [""] * (width - len(r)) for r in rows]
    return pad(head), pad(body)


def header_labels(head_rows, width):

    if not head_rows:
        return [""] * width
    labels = []
    for c in range(width):
        parts = []
        for row in head_rows:
            if c < len(row):
                v = row[c].strip()
                if v and (not parts or parts[-1] != v):
                    parts.append(v)
        labels.append(" | ".join(parts))
    return labels


def iter_table_wraps(root):

    for tw in root.iter("table-wrap"):
        label = ""
        lab = tw.find("label")
        if lab is not None:
            label = cell_text(lab)
        caption = ""
        cap = tw.find("caption")
        if cap is not None:
            caption = cell_text(cap)
        footer = ""
        wrapfoot = tw.find("table-wrap-foot")
        if wrapfoot is not None:
            footer = cell_text(wrapfoot)

        table = tw.find(".//table")
        if table is None:
            continue
        head, body = build_grid(table)
        width = max([len(r) for r in head + body] or [0])
        yield {
            "id": tw.get("id") or "",
            "label": label,
            "caption": caption,
            "footer": footer,
            "head": head,
            "body": body,
            "width": width,
            "headers": header_labels(head, width),
        }
