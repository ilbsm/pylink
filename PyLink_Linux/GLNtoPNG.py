#!/usr/bin/python
import sys
import math
from PIL import Image, ImageDraw, ImageFont

EPS = 0.0001
g = 0.69  # TRESHOLDS for lasso classification
h = 0.55
sl = 1.5
R = 2  # precision for rounding an output

WIDTH = 400
HEIGHT = 400
FONT_SIZE = WIDTH / 27
FONT_SIZE_SMALL = WIDTH / 35
FONT_TYPE = "/usr/share/fonts/truetype/freefont/FreeMono.ttf"
FONT_TYPE_SMALL = "/usr/share/fonts/truetype/freefont/FreeMono.ttf"


def compareEq(a, b):
    if not type(a) == list and not type(a) == tuple:
        return abs(a - b) < EPS
    else:
        assert len(a) == len(b)
        for i in range(len(a)):
            if not compareEq(a[i], b[i]): return False
        return True


def compareGeq(a, b): return a > b - EPS


def compareGt(a, b): return a > b + EPS


def scalarProduct(a, b):
    # return scalar product of two vectors (assumption: a,b vectors (lists) of three double numbers
    assert (type(a) == list and type(b) == list and len(a) == 3 and len(b) == 3)
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])


def vectorProduct(a, b):
    # return vector product of two vectors (assumption: a,b vectors (lists) of three double numbers
    assert (type(a) == list and type(b) == list and len(a) == 3 and len(b) == 3)
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def normalizedVectorProduct(a, b):
    # return vp := a x b / |a x b| (assumption: a,b vectors (lists) of three double numbers
    assert (type(a) == list and type(b) == list and len(a) == 3 and len(b) == 3)
    v = vectorProduct(a, b)
    d = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** 0.5
    if d == 0: return 0
    return [el / d for el in v]


def fileCheck(fn):
    try:
        open(fn, "r")
        return 1
    except IOError:
        print(("Error: File", fn, "does not appear to exist.\n"))
        return 0


def chainRead(filename, begin=False, end=False):
    # Read coordinates from .xyz file, cuts them from id "begin" to "end", and put them into the returned list; empty list -> error/problem -> end of the program; works with 4 and 5 columns file.
    if (not fileCheck(filename)): return []
    comp = []
    b, e = begin, end
    if begin != False and begin > end: b, e = end, begin

    k, x, y, z = -1, -1, -1, -1
    f = open(filename)
    for line in f.readlines():
        if len(line) > 1:
            k2, x2, y2, z2 = [float(x) for x in line.split()[:4]]
            if compareEq([k, x, y, z], [k2, x2, y2, z2]):
                print(("Error: There were two identical atoms with the index", k,
                      "in the file. Program is done becouse of that.\n"))
                return []
            k, x, y, z = k2, x2, y2, z2
            if (k >= b and k <= e) or (not begin and not end):
                comp.append((int(k), (x, y, z)))
    f.close()
    return comp


def CompsOverlie(comp1, comp2):
    # Return True if there are points on two comps that overlie
    for i in range(len(comp1)):
        for j in range(len(comp2)):
            if compareEq(comp1[i][1], comp2[j][1]): return True
    return False


def linking_oneSegment(A1, A2, B1, B2):
    # returns self linking number between A1A2 and B1B2
    # Ken's notations

    a = [B1[i] - A1[i] for i in [0, 1, 2]]
    b = [B2[i] - A1[i] for i in [0, 1, 2]]
    c = [B2[i] - A2[i] for i in [0, 1, 2]]
    d = [B1[i] - A2[i] for i in [0, 1, 2]]

    n1 = normalizedVectorProduct(a, b)
    n2 = normalizedVectorProduct(b, c)
    n3 = normalizedVectorProduct(c, d)
    n4 = normalizedVectorProduct(d, a)
    if n1 == 0 or n2 == 0 or n3 == 0 or n4 == 0:
        print("Problem with Vector Product...", A1, A2, B1, B2, a, b, c, d)
        return 0

    f = [B2[i] - B1[i] for i in [0, 1, 2]]
    h = [A2[i] - A1[i] for i in [0, 1, 2]]
    s = vectorProduct(f, h)

    x = scalarProduct(s, a)
    sign = 0
    if compareEq(x, 0): x, sign = 0, 0
    if x > 0: sign = 1
    if x < 0: sign = -1

    as1, as2, as3, as4 = scalarProduct(n1, n2), scalarProduct(n2, n3), scalarProduct(n3, n4), scalarProduct(n4, n1)
    if (compareGeq(as1, 1)): as1 = 1
    if (compareGeq(as2, 1)): as2 = 1
    if (compareGeq(as3, 1)): as3 = 1
    if (compareGeq(as4, 1)): as4 = 1
    if (compareGeq(-1, as1)): as1 = -1
    if (compareGeq(-1, as2)): as2 = -1
    if (compareGeq(-1, as3)): as3 = -1
    if (compareGeq(-1, as4)): as4 = -1

    return sign * (math.asin(as1) + math.asin(as2) + math.asin(as3) + math.asin(as4)) / (4 * math.pi)


def linking_components(comp1, comp2):
    N1, N2 = len(comp1), len(comp2)
    links = []  # links[i][j] = GLN ( comp1<i,i+1>, comp2<j,j+1> ) - GLN between each pair of segments: ONE SEGMENT with ONE SEGMENT
    #    links = (N1-1)*[(N2-1)*[0]]
    for i in range(N1 - 1): links.append((N2 - 1) * [0])

    for i in range(N1 - 1):
        for j in range(N2 - 1):
            links[i][j] = linking_oneSegment(comp1[i][1], comp1[i + 1][1], comp2[j][1], comp2[j + 1][1])

    vlinksLOOP1, vlinksLOOP2 = (N2 - 1) * [0], (N1 - 1) * [
        0]  # vlinksLOOP1[j] = GLN ( comp1, comp2<j,j+1> ) - GLN between WHOLE comp1 (as a LOOP1) and ONE SEGMENT of comp2
    for i in range(N1 - 1):
        for j in range(N2 - 1):
            vlinksLOOP1[j] += links[i][j]
            vlinksLOOP2[i] += links[i][j]

    linksLOOP1, linksLOOP2 = [], []  # linksLOOP1[i][j] = GLN ( comp1, comp2<i-j> ) - GLN between WHOLE comp1 (as a LOOP1) and PART of comp2 <i-j>
    for j in range(N2): linksLOOP1.append(N2 * [0])
    for i in range(N1): linksLOOP2.append(N1 * [0])

    for j1 in range(N2):
        for j2 in range(j1 + 1, N2):
            linksLOOP1[j1][j2] = linksLOOP1[j1][j2 - 1] + vlinksLOOP1[j2 - 1]
    for i1 in range(N1):
        for i2 in range(i1 + 1, N1):
            linksLOOP2[i1][i2] = linksLOOP2[i1][i2 - 1] + vlinksLOOP2[i2 - 1]

    return linksLOOP1, linksLOOP2


def colorFromGLN(gln):
    if (gln < -1):
        return (int(255 * 1 / (gln * gln)), 0, 0)
    elif (gln <= 0):
        return (255, int(255 * (1 + gln)), int(255 * (1 + gln)))
    elif (gln <= 1):
        return (int(255 * (1 - gln)), int(255 * (1 - gln)), 255)
    else:
        return (0, 0, int(255 * 1 / (gln * gln)))


def writePng(linksLOOP, filename, xxx_todo_changeme, xxx_todo_changeme1):
    (mmax, xmax, ymax) = xxx_todo_changeme
    (mmin, xmin, ymin) = xxx_todo_changeme1
    CHAIN = len(linksLOOP)
    KLATKA = int(WIDTH / CHAIN)
    FLOATKLATKA = float(WIDTH) / CHAIN

    image = Image.new('RGB', (WIDTH, HEIGHT), 'white')
    draw = ImageDraw.Draw(image)

    for i in range(CHAIN):
        for j in range(i + 1, CHAIN):
            color = colorFromGLN(linksLOOP[i][j])
            for x in range(0, KLATKA + 1):
                for y in range(0, KLATKA + 1):
                    draw.point((int(i * FLOATKLATKA) + x, int(j * FLOATKLATKA) + y), fill=color)

    # Extrema
    font = ImageFont.truetype(FONT_TYPE, FONT_SIZE)
    color_min, color_max = colorFromGLN(linksLOOP[xmin][ymin]), colorFromGLN(linksLOOP[xmax][ymax])

    xmax = int(xmax * FLOATKLATKA) if int(xmax * FLOATKLATKA) < HEIGHT - FONT_SIZE else HEIGHT - FONT_SIZE - 5
    ymax = int(ymax * FLOATKLATKA) if int(ymax * FLOATKLATKA) < WIDTH - FONT_SIZE else WIDTH - FONT_SIZE - 5
    xmin = int(xmin * FLOATKLATKA) if int(xmin * FLOATKLATKA) < HEIGHT - FONT_SIZE else HEIGHT - FONT_SIZE - 5
    ymin = int(ymin * FLOATKLATKA) if int(ymin * FLOATKLATKA) < WIDTH - FONT_SIZE else WIDTH - FONT_SIZE - 5

    r = int(0.4 * FLOATKLATKA)
    if CHAIN > 25: r = int(0.6 * FLOATKLATKA)
    if CHAIN > 50: r = int(0.9 * FLOATKLATKA)
    if CHAIN > 80: r = int(1.5 * FLOATKLATKA)
    if CHAIN > 120: r = int(2 * FLOATKLATKA)
    if CHAIN > 200: r = int(4 * FLOATKLATKA)
    if CHAIN > 500: r = int(8 * FLOATKLATKA)
    draw.ellipse((xmin, ymin, xmin + 2 * r, ymin + 2 * r), fill=color_min, outline=(0, 0, 0, 255))
    draw.ellipse((xmax, ymax, xmax + 2 * r, ymax + 2 * r), fill=color_max, outline=(0, 0, 0, 255))
    draw.ellipse((0.550 * WIDTH, 0.09 * WIDTH, 0.590 * WIDTH, 0.13 * WIDTH), fill=color_min, outline=(0, 0, 0, 255))
    draw.ellipse((0.550 * WIDTH, 0.15 * WIDTH, 0.590 * WIDTH, 0.19 * WIDTH), fill=color_max, outline=(0, 0, 0, 255))

    draw.text((0.60 * WIDTH, 0.09 * WIDTH), "min GLN = " + str(mmin), font=font, fill=(0, 0, 0, 255))
    draw.text((0.60 * WIDTH, 0.15 * WIDTH), "max GLN = " + str(mmax), font=font, fill=(0, 0, 0, 255))

    # Osie
    font = ImageFont.truetype(FONT_TYPE_SMALL, FONT_SIZE_SMALL)
    step = int(CHAIN / 5) / 10 * 10 if CHAIN > 150 else int(CHAIN / 5) / 5 * 5 if CHAIN > 50 else int(CHAIN / 5)
    for x in range(1, 6):
        draw.text((x * step * FLOATKLATKA, HEIGHT - FONT_SIZE_SMALL - 2), str(x * step), font=font,
                  fill=(50, 50, 50, 255))
        draw.text((0, x * step * FLOATKLATKA), "-" + str(x * step), font=font, fill=(50, 50, 50, 255))

    image.save(filename + '.png')


def main():
    # Reading arguments
    if (len(sys.argv) == 1):
        print(("*\nUsage of the program: python " + sys.argv[0] + " <filename1> (<begin of a comp1> <end of a comp1>) "
                                                                 "<filename2> (<begin of a comp2> <end of a comp2>) "
                                                                 "(<-additional_option and argument>^n)"))
        print("*\nAdditional options:\n   -> -file (0,1): if we create png file/s; implicitly file=1;\n   "
              "-> -out (0,1,2): output on the screen, 0:none, 1:short, 2:long<compatybile with Wanda's GLN paper>; "
              "implicitly out=1;\n   -> -close (00,10,01,11): if we close 1st/2nd comp; implicitly close=00;\n   "
              "-> -loop (0,1,2): if we calculate bots \"glns\", 0:both, 1:only 1st comp as a loop, 2:only 2nd comp "
              "as a loop; implicitly loop=0;\n   -> -size (positive int<10000): size of the png picture; implicitly "
              "size=WIDTH(at top)=500 for now;\n   -> -ht (string): if you want to have that string at the beginning "
              "of png files.\n*")
        return 0

    b1, e1, b2, e2 = False, False, False, False
    ar = 1
    filename1 = sys.argv[ar]
    ar = ar + 1
    if (sys.argv[ar]).isdigit():
        b1, e1 = int(sys.argv[ar]), int(sys.argv[ar + 1])
        ar = ar + 2
    filename2 = sys.argv[ar]
    ar = ar + 1
    if len(sys.argv) > ar and (sys.argv[ar]).isdigit():
        b2, e2 = int(sys.argv[ar]), int(sys.argv[ar + 1])
        ar = ar + 2

    # Read optional arguments
    pngfile = True  # if we create pngfile/s
    out = 1  # output on the screen, 0-none, 1-short, 2-long
    cl1, cl2 = False, False  # if we threat 1/2 comp as closed
    loop = 0  # if we check both combinations - 0-both, 1-comp1 as a loop, 2-comp2 as a loop
    hashtag = ""

    while ar < len(sys.argv) - 1:  # -file, -out, -closed, -loop
        if sys.argv[ar] == "-file":
            if sys.argv[ar + 1] == "0": pngfile = False
            ar += 2
        elif sys.argv[ar] == "-out":
            out = int(sys.argv[ar + 1])
            if out < 0 or out > 2: out = 1
            ar += 2
        elif sys.argv[ar] == "-close":
            if len(sys.argv[ar + 1]) != 2:
                print("Error: argument after -close must consists of two digits (0-1).")
                return 0
            if sys.argv[ar + 1][0] == "1": cl1 = True
            if sys.argv[ar + 1][1] == "1": cl2 = True
            ar += 2
        elif sys.argv[ar] == "-loop":
            loop = int(sys.argv[ar + 1])
            if loop < 0 or loop > 2: loop = 0
            ar += 2
        elif sys.argv[ar] == "-size":
            size = int(sys.argv[ar + 1])
            if size < 10000 or size > 50: WIDTH, HIGHT = size, size
            ar += 2
        elif sys.argv[ar] == "-ht":
            hashtag = sys.argv[ar + 1] + "_"
            ar += 2
        elif sys.argv[ar] == "-ion":
            if int(sys.argv[ar + 1]) == 1: b1, e1, b2, e2 = False, False, False, False
            ar += 2

    # Creating components
    comp1 = chainRead(filename1, b1, e1)
    if cl1: comp1.append(comp1[0])
    comp2 = chainRead(filename2, b2, e2)
    if cl2: comp2.append(comp2[0])
    N1, N2 = len(comp1), len(comp2)

    if N1 == 0 or N2 == 0:
        print("empty-comp . . . . . . . . . . . . . . .")
        return 0
    if N1 == 1 or N2 == 1:
        print("1-length-comp . . . . . . . . . . . . . . .")
        return 0
    if CompsOverlie(comp1, comp2):
        print("comps-overlie . . . . . . . . . . . . . . .")
        return 0

    # Calculate linking
    linksLOOP1, linksLOOP2 = linking_components(comp1, comp2)

    # Preparing data for out
    #    out =1: wh1 max1 min1 cl1
    #    out =2: spacja N2 wh1: wh1+: max1: min1: max10: min10: cl

    whole1, whole2 = linksLOOP1[0][N2 - 1], linksLOOP2[0][N1 - 1]
    max1, xmax1, ymax1, min1, xmin1, ymin1, max2, xmax2, ymax2, min2, xmin2, ymin2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    max10, ymax10, min10, ymin10, max20, ymax20, min20, ymin20 = 0, 0, 0, 0, 0, 0, 0, 0

    if loop == 0 or loop == 2:
        for i1 in range(N1):
            for i2 in range(i1 + 1, N1):
                if linksLOOP2[i1][i2] > max2: max2, xmax2, ymax2 = linksLOOP2[i1][i2], i1, i2
                if linksLOOP2[i1][i2] < min2: min2, xmin2, ymin2 = linksLOOP2[i1][i2], i1, i2
        if linksLOOP2[0][i1] > max20: max20, ymax20 = linksLOOP2[0][i1], i1
        if linksLOOP2[0][i1] < min20: min20, ymin20 = linksLOOP2[0][i1], i1
    if loop == 0 or loop == 1:
        for j1 in range(N2):
            for j2 in range(j1 + 1, N2):
                if linksLOOP1[j1][j2] > max1: max1, xmax1, ymax1 = linksLOOP1[j1][j2], j1, j2
                if linksLOOP1[j1][j2] < min1: min1, xmin1, ymin1 = linksLOOP1[j1][j2], j1, j2
        if linksLOOP1[0][j1] > max10: max10, ymax10 = linksLOOP1[0][j1], j1
        if linksLOOP1[0][j1] < min10: min10, ymin10 = linksLOOP1[0][j1], j1

    whole1, whole2 = round(whole1, R), round(whole2, R)
    max1, min1, max2, min2, max10, min10, max20, min20 = round(max1, R), round(min1, R), round(max2, R), round(min2,
                                                                                                               R), round(
        max10, R), round(min10, R), round(max20, R), round(min20, R)

    amin1, amin2, awhole1, awhole2 = abs(min1), abs(min2), abs(whole1), abs(whole2)
    if max(max1, amin1) >= sl:
        cl1 = "LS"
    elif max(max1, amin1) < g:
        cl1 = "L0"
    elif amin1 < g <= max1 or amin1 >= g > max1:
        cl1 = "L1"
    elif awhole1 < h:
        cl1 = "L2+"
    else:
        cl1 = "L3+"
    if max(max2, amin2) >= sl:
        cl2 = "LS"
    elif max(max2, amin2) < g:
        cl2 = "L0"
    elif amin2 < g <= max2 or amin2 >= g > max2:
        cl2 = "L1"
    elif awhole2 < h:
        cl2 = "L2+"
    else:
        cl2 = "L3+"

    if out == 1:
        if loop == 0 or loop == 1:
            print(( "*LOOP1 wh: " + str(whole1) + " max: " + str(max1) + " min: " + str(min1) + " cl: " + cl1))
        if loop == 0 or loop == 2:
            print(( "*LOOP2 wh: " + str(whole2) + " max: " + str(max2) + " min: " + str(min2) + " cl: " + cl2))
    if out == 2:
        if loop == 0 or loop == 1:
            print((" " + str(N2) + " wh1: " + str(whole1) + " wh1+: 0 max1: " + str(max1) + " min1: " + str(min1) +
                  " max10: " + str(max10) + " " + str(ymax10) + " min10: " + str(min10) + " " + str(ymin10) +
                  " cl1: " + cl1))
        if loop == 0 or loop == 2:
            print((" " + str(N1) + " wh2: " + str(whole2) + " wh2+: 0 max2: " + str(max2) + " min2: " + str(min2) +
                  " max20: " + str(max20) + " " + str(ymax20) + " min20: " + str(min20) + " " + str(ymin20) +
                  " cl2: " + cl2))

    i1, i2 = filename1.rfind("/"), filename1.rfind(".")
    n1 = filename1[i1 + 1:] if i2 == -1 else filename1[
                                             i1 + 1:i2]  # n1 - just the name of the first file, without a path and an extension
    if b1 != False: n1 += "_L" + str(b1) + "-" + str(e1)
    j1, j2 = filename2.rfind("/"), filename2.rfind(".")
    n2 = filename2[j1 + 1:] if j2 == -1 else filename2[
                                             j1 + 1:j2]  # n2 - just the name of the second file, without a path and an extension
    if b2 != False: n2 += "_L" + str(b2) + "-" + str(e2)

    name1 = "" if i1 == -1 else filename1[:i1 + 1]  # paths
    name2 = "" if j1 == -1 else filename2[:j1 + 1]
    if "macro_" in n1 and "macro_" in n2:
        n1 = "1" + n1.split("macro_")[1].split("_L")[0]
        n2 = "2" + n2.split("macro_")[1].split("_L")[0]
    else:
        n1 = "1" + n1.split(".pdb_")[1][0]
        n2 = "2" + n2.split(".pdb_")[1][0]
    name1 = name1 + hashtag + "GLN_" + n1 + "-" + n2
    name2 = name2 + hashtag + "GLN_" + n2 + "-" + n1

    if pngfile:
        if loop == 0 or loop == 1: writePng(linksLOOP1, name1, (max1, xmax1, ymax1), (min1, xmin1, ymin1))
        if loop == 0 or loop == 2: writePng(linksLOOP2, name2, (max2, xmax2, ymax2), (min2, xmin2, ymin2))


##########################
if __name__ == "__main__":
    main()
