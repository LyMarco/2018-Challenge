import svgwrite
from setUpCh import chSetUp
from prepTADData import prepData
from prepGeneData import prepGeneData
from math import cos, sin, pi

CENTER = (500, 150)      #where to center our circos plot
RADIUS = 125            #size of circos plot
SUBCENTER = (500, 500)
SUBRADIUS = 150
START = (0, RADIUS)     #for calculating gene arcs
SUBSTART = (0, SUBRADIUS)
CHROMLEN = 64444167     #nts on chromosome 20

#gene location and GOA based print - circos
#draw the chromosome as a circle and add arcs to represent gene occupancy
#add edges between genes that share the same GOA

def calcArcCoords(chDF, starting, centered):
    #find the coordinates of the start of the arc
    chDF["thetaStart"] = (chDF.start/chDF.circle_len)*(2*pi)
    chDF["arcStartX"] = chDF.thetaStart.apply(cos) * starting[0] - chDF.thetaStart.apply(sin) * starting[1]
    chDF["arcStartY"] = chDF.thetaStart.apply(sin) * starting[0] + chDF.thetaStart.apply(cos) * starting[1]

    #find the coordinates of the end of the arc
    chDF["thetaEnd"] = (chDF.end/chDF.circle_len)*(2*pi)
    chDF["arcEndX"] = chDF.thetaEnd.apply(cos) * starting[0] - chDF.thetaEnd.apply(sin) * starting[1]
    chDF["arcEndY"] = chDF.thetaEnd.apply(sin) * starting[0] + chDF.thetaEnd.apply(cos) * starting[1]

    #find the coordinates of the middle of the arc(for edges to originate from)
    chDF["thetaMid"] = ((chDF.start+((chDF.end-chDF.start)/2))/chDF.circle_len)*(2*pi)
    chDF["arcMidX"] = chDF.thetaMid.apply(cos) * starting[0] - chDF.thetaMid.apply(sin) * starting[1]
    chDF["arcMidY"] = chDF.thetaMid.apply(sin) * starting[0] + chDF.thetaMid.apply(cos) * starting[1]

    #now shift all coordinates to the defined center
    chDF.arcStartX = chDF.arcStartX + centered[0]
    chDF.arcStartY = chDF.arcStartY + centered[1]
    chDF.arcEndX = chDF.arcEndX + centered[0]
    chDF.arcEndY = chDF.arcEndY + centered[1]
    chDF.arcMidX = chDF.arcMidX + centered[0]
    chDF.arcMidY = chDF.arcMidY + centered[1]

    return chDF


def circoGeneDraw(chDF, genesDF, fName):

    dwg = svgwrite.Drawing(filename=fName)

    #add title
    dwg.add(dwg.text("CHR 20 TADS", insert=CENTER, fill="black", text_anchor="middle"))

    #draw chromosome line
    dwg.add(dwg.circle(CENTER, RADIUS, fill_opacity=0.0, stroke="black", stroke_width=1))

    from_ = "0 " + str(CENTER[0]) + " " + str(CENTER[1])
    to_ = "360 " + str(CENTER[0]) + " " + str(CENTER[1])

    #draw arcs for each tad
    tadArcs = dwg.add(dwg.g(id="tadArcs", fill="white", stroke_width=5))
    tadArcs.add(
        dwg.animateTransform("rotate", "transform", id="arcs", from_=from_, to=to_, dur="360s",
                             begin="0s", repeatCount="indefinite"))
    for g in chDF.index:
        gStart = "M" + str(chDF.loc[g].arcStartX) + "," + str(chDF.loc[g].arcStartY)
        gArc = "A" + str(RADIUS) + "," + str(RADIUS)
        gEnd = str(chDF.loc[g].arcEndX) + "," + str(chDF.loc[g].arcEndY)
        tadArcs.add(dwg.path(d=gStart + gArc + " 0 0,1 " + gEnd, stroke=chDF.loc[g].colour))

        end_angle = 360 - chDF.loc[g].start/CHROMLEN * 360
        start_angle = 360 - chDF.loc[g].end/CHROMLEN * 360

        start_time = str(start_angle) + "s"
        end_time = str(end_angle) + "s"

        sub_circle = dwg.add(dwg.circle(SUBCENTER, SUBRADIUS, opacity=0, stroke="black", stroke_width=1,
                                        fill=chDF.loc[g].colour))
        sub_circle.add(dwg.animate(attributeName="opacity", id="fills", from_="1", to="0",
                                                begin=start_time, end=end_time, dur="360s", repeatCount="indefinite"))

        geneArcs = dwg.add(dwg.g(id="geneArcs", fill="white", stroke_width=5))


        sub_genes_df = genesDF[genesDF["TAD"] == g]

        for h in sub_genes_df.index:

            # look up the gene arc start and end points and plug into drawing
            hStart = "M" + str(sub_genes_df.loc[h].arcStartX) + "," + str(sub_genes_df.loc[h].arcStartY)
            hArc = "A" + str(SUBRADIUS) + "," + str(SUBRADIUS)
            hEnd = str(sub_genes_df.loc[h].arcEndX) + "," + str(sub_genes_df.loc[h].arcEndY)
            geneArcs.add(dwg.path(d=hStart + hArc + " 0 0,1 " + hEnd, opacity=0, stroke=chDF.loc[g].colour))
            geneArcs.add(dwg.animate(attributeName="opacity", id="fill", from_="1", to="0",
                                       begin=start_time, end=end_time, dur="360s", repeatCount="indefinite"))
        # for h in genesDF.index:
        #     if genesDF.loc[h].TAD == g:
        #         # look up the gene arc start and end points and plug into drawing
        #         hStart = "M" + str(genesDF.loc[h].arcStartX) + "," + str(genesDF.loc[h].arcStartY)
        #         hArc = "A" + str(SUBRADIUS) + "," + str(SUBRADIUS)
        #         hEnd = str(genesDF.loc[h].arcEndX) + "," + str(genesDF.loc[h].arcEndY)
        #         geneArcs.add(dwg.path(d=hStart + hArc + " 0 0,1 " + hEnd, stroke=chDF.loc[g].colour, opacity=0))
        #         geneArcs.add(dwg.animate(attributeName="opacity", id="fills", from_="1", to="0",
        #                                         begin=start_time, end=end_time, dur="360s", repeatCount="indefinite"))
    dwg.save()


if __name__ == '__main__':
    print("setting up...")
    ch20 = prepData()
    geneData = prepGeneData()
    print("set up complete\ncomputing gene arcs...")
    ch20["circle_len"] = CHROMLEN
    ch20 = calcArcCoords(ch20, START, CENTER)
    geneData = calcArcCoords(geneData, SUBSTART, SUBCENTER)
    # print("gene arcs computed\ncomputing edges... (this may take a few seconds)")
    # ch20Edges = makeEdgeList(ch20)
    print("edges computed\nrendering svg... (this may take a few seconds)")
    circoGeneDraw(ch20, geneData, "circosDraw.svg")

