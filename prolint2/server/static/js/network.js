import { radarApp } from "../js/radialApp.js";
// ##########################
export function networkApp(lipid = "CHOL") {

    var root = am5.Root.new("chartdiv_y");

    root.setThemes([
        am5themes_Animated.new(root)
    ]);

    var linkDefaultOpacity = 0.5,
        lihnkHoveredOpacity = 1;

    // Create series
    var series = root.container.children.push(am5flow.ChordNonRibbon.new(root, {
        sourceIdField: "from",
        targetIdField: "to",
        valueField: "value",
        padAngle: 0,
        startAngle: 90,
        draggable: true,
        lipid: lipid,

    }));

    series.nodes.labels.template.setAll({
        textType: "radial",
        fontSize: "0.6em",
        radius: 15,
    });

    series.nodes.bullets.push(function (_root, _series, dataItem) {
        var bulletCircle = am5.Circle.new(root, {
            radius: 3.5,
            fill: am5.color("#8E8A8A"),
            fillOpacity: 1,
        });

        bulletCircle.adapters.add("fill", function (fill, target) {
            var dataItem = target.dataItem;
            if (dataItem) {
                var sum = dataItem.get("sum", 0) // we can also use sumIncoming or sumOutgoing
                var min = Infinity;
                var max = -Infinity;
                am5.array.each(series.nodes.dataItems, function (dataItem) {
                    var value = dataItem.get("sum");
                    if (value < min) {
                        min = value;
                    }
                    if (value > max) {
                        max = value;
                    }
                })

                return am5.Color.interpolate((sum - min) / (max - min), am5.color(0x673AB7), am5.color(0xF44336));
            }
            return fill;
        })

        return am5.Bullet.new(root, {
            sprite: bulletCircle
        });
    });

    series.children.moveValue(series.bulletsContainer, 0);

    var hoverColor = am5.color(0x393D47); 
    // am5.array.each(subSeries.dataItems, function (dataItem, ix) {
    //     if (dataItem.dataContext.category == lipid) {
    //         hoverColor = subSeries.get("colors").getIndex(ix)
    //     }
    // })

    // Add heat legend
    var networkHeatLegend = series.children.push(am5.HeatLegend.new(root, {
        orientation: "horizontal",
        startColor: am5.color(0x673AB7),
        endColor: am5.color(0xF44336),
        width: 160,
        x: am5.percent(70),
        y: am5.percent(95),
        opacity: 1,
    }));

    series.nodes.rectangles.template.events.on("pointerover", function (ev) {
        var di = ev.target.dataItem;
        if (di) {
            networkHeatLegend.showValue(di._settings.sum);
        }

        var incomingLinks = di.get('incomingLinks')
        var outgoingLinks = di.get('outgoingLinks')
        var selectSections = [{
            residue_number: parseInt(di.dataContext.id),
            color: {
                r: 255,
                g: 0,
                b: 255
            }
        }]
        if (incomingLinks) {
            incomingLinks.forEach(link => {
                var sourceId = link.get('sourceId')
                if (sourceId) {
                    selectSections.push({
                        residue_number: parseInt(sourceId),
                        color: {
                            r: 255,
                            g: 0,
                            b: 255
                        }
                    })
                }
                link._settings.link._settings.stroke = hoverColor
                link._settings.link._display.alpha = lihnkHoveredOpacity
            })
        }
        if (outgoingLinks) {
            outgoingLinks.forEach(link => {
                var targetId = link.get('targetId')
                selectSections.push({
                    residue_number: parseInt(targetId),
                    color: {
                        r: 255,
                        g: 0,
                        b: 255
                    }
                })
                link._settings.link._settings.stroke = hoverColor
                link._settings.link._display.alpha = lihnkHoveredOpacity
            })
        }
        viewerInstance.visual.highlight({
            data: selectSections,
        })
    })

    series.nodes.rectangles.template.events.on("pointerout", function (ev) {

        var incomingLinks = ev.target._dataItem._settings.incomingLinks
        var outgoingLinks = ev.target.dataItem._settings.outgoingLinks
        if (incomingLinks) {
            incomingLinks.forEach(link => {
                link._settings.link._settings.stroke = am5.color("#8E8A8A")
                link._settings.link._display.alpha = linkDefaultOpacity
            })
        }
        if (outgoingLinks) {
            outgoingLinks.forEach(link => {
                link._settings.link._settings.stroke = am5.color("#8E8A8A")
                link._settings.link._display.alpha = linkDefaultOpacity
            })
        }
        var selectSections = []
        viewerInstance.visual.highlight({
            data: selectSections,
        })
    })

    series.nodes.rectangles.template.events.on('click', function (ev) {

        var di = ev.target.dataItem;
        var incomingLinks = di.get('incomingLinks')
        var outgoingLinks = di.get('outgoingLinks')
        var selectSections = [{
            residue_number: parseInt(di.dataContext.id),
            representationColor: {r:255, g:0, b:0},
            representation: 'spacefill',
        }]
        if (incomingLinks) {
            incomingLinks.forEach(link => {
                var sourceId = link.get('sourceId')
                if (sourceId) {
                    selectSections.push({
                        residue_number: parseInt(sourceId),
                        representationColor: {r:255, g:255, b:255},
                        representation: 'spacefill',
                    })
                }
            })
        }
        if (outgoingLinks) {
            outgoingLinks.forEach(link => {
                var targetId = link.get('targetId')
                selectSections.push({
                    residue_number: parseInt(targetId),
                    representation: 'spacefill',
                })
            })
        }
        viewerInstance.visual.select({
            data: selectSections,
        })

    })

    obj = {
        "lipid": lipid,
        // "residueID": ctx.category
    }
    fetch('/network/' + JSON.stringify(obj))
        .then(response => response.json())
        .then(responseData => {

            var data = responseData['chordElements']
            var posRes = responseData['positionResidues']
            
            series.data.setAll(data);
            series.lipidNodes = responseData['lipidNodes'];

            series.links.template.setAll({
                strokeWidth: 0.2,
                opacity: linkDefaultOpacity,
                stroke: am5.color("#8E8A8A"),
                strokeStyle: "none",
            })

            series.nodes.rectangles.template.setAll({
                fillOpacity: 0,
                fill: hoverColor,
                tooltipText: "Residue [bold]{name}[/]\nShared Contacts: {sum}"
            });

            series.events.on("datavalidated", function (ev) {
                var el = ev.target;

                var bulletSums = [];
                for (let ix = 0; ix < series.nodes.dataItems.length; ix++) {
                    var di = series.nodes.dataItems[ix];
                    var sum = di.get('sum');
                    if (sum != 0) {
                        bulletSums.push(sum);
                    }
                }

                networkHeatLegend.set("startValue", Math.min(...bulletSums));
                networkHeatLegend.set("endValue", Math.max(...bulletSums));

                for (let jx = 0; jx < posRes.length; jx++) {
                    const pos = posRes[jx];
                    el.nodes.dataItems[pos].bullets[0]._settings.sprite._display.visible = false
                    el.nodes.labels._values[pos]._display.visible = false
                    el.nodes.rectangles.template._entities[pos]._display.visible = false
                }
                el.links._values.forEach((nodeLinks, ix) => {
                    if (nodeLinks.dataItem.dataContext.from == 0) {
                        nodeLinks._settings.strokeWidth = 0;
                        nodeLinks._settings.strokeOpacity = 0;
                        nodeLinks._display.visible = false
                    } else {
                        nodeLinks._settings.strokeWidth = nodeLinks.dataItem.dataContext.valueWidth * 2;
                    }
                })
            });
        });

    series.appear(1000, 100);

    return {
        "series": series,
        "root": root,
        "active": true
    }
}

// export default networkApp