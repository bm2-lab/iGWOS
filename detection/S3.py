import svgwrite


def SITE_visualizeOfftargets(offtargets, ref_seq, outfile, title=None,box_size=15,v_spacing=3,colors='default'):
    if colors =='default':
        colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', '-': '#B3B3B3'}

    dwg = svgwrite.Drawing(outfile + '.svg', profile='full', size=(u'100%', 100 + len(offtargets)*(box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # Draw ticks
    tick_locations = [1, len(ref_seq)]
    tick_locations += range(len(ref_seq) + 1)[::10][1:]
    for x in tick_locations:
        dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(ref_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 10  # leave some extra space after the reference row
    for j, seq in enumerate(offtargets):
        realigned_ref_seq =  offtargets[j]['realigned_ref_seq']
        k = 0
        y = y_offset + j * box_size
        for i, (c, r) in enumerate(zip(seq['seq'], realigned_ref_seq)):
            x = x_offset + k * box_size
            if r == '-':
                if 0 < k < len(realigned_ref_seq):
                    x = x_offset + (k - 0.25) * box_size
                    dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))

            elif c == r or r == 'N':
                 dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                 k += 1
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                k += 1


        reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text)

    dwg.save()