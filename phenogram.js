function renderTaxon (taxon, tree)
{
    if (taxon["subtaxa"])
        for (subtaxon of taxon.subtaxa)
        {
            let line = document.createElementNS("http://www.w3.org/2000/svg", "line");
            line.setAttribute("x1", 100 * taxon.depth + "%");
            line.setAttribute("y1", 100 * taxon.index + "%");
            line.setAttribute("x2", 100 * subtaxon.depth + "%");
            line.setAttribute("y2", 100 * subtaxon.index + "%");
            line.setAttribute("stroke", "black");
            line.setAttribute("stroke-width", "2");
            tree.appendChild(line);
            for (const subtaxon of taxon.subtaxa)
                renderTaxon(subtaxon, tree);
        }
    else
    {
        line = document.createElementNS("http://www.w3.org/2000/svg", "line");
        line.setAttribute("x1", 100 * taxon.depth + "%");
        line.setAttribute("y1", 100 * taxon.index + "%");
        line.setAttribute("x2", 100 * 1 + "%");
        line.setAttribute("y2", 100 * taxon.index + "%");
        line.setAttribute("stroke", "gray");
        tree.appendChild(line);
    }
    
    let circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
    circle.setAttribute("cx", 100 * taxon.depth + "%");
    circle.setAttribute("cy", 100 * taxon.index + "%");
    circle.setAttribute("r", 3);
    tree.appendChild(circle);
}

function renderExpression(taxon, element)
{
    if (taxon["subtaxa"])
        for (subtaxon of taxon.subtaxa)
            renderExpression(subtaxon, element);
    else
    {
        let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
        svg.setAttribute("height", 20);
        for (let i = 0; i < taxon.markers.length; ++i)
        {
            let circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
            circle.setAttribute("cx", 10 + i * 20);
            circle.setAttribute("cy", 10);
            circle.setAttribute("r", 5);
            svg.appendChild(circle);
    
        }
        element.appendChild(svg);
        element.appendChild(document.createElement("br"));
    }
}


function render (taxonomy)
{
    let svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    renderTaxon(taxonomy, svg);

    let div = document.createElement("div");
    renderExpression(taxonomy, div);

    document.getElementById("taxonomy_tree").appendChild(svg);
    document.getElementById("taxonomy_expression").appendChild(div);
}
