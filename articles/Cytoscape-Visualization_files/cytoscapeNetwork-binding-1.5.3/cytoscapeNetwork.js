/* ==========================================================================
   cytoscapeNetwork.js
   htmlwidgets binding for the cytoscapeNetwork R package.

   Dependencies (declared in cytoscapeNetwork.yaml):
     - cytoscape.min.js
     - dagre.min.js / graphlib.min.js
     - cytoscape-dagre.js

   The `x` object passed from R (via createWidget / jsonlite serialisation):
   {
     nodes        : [ { id, label, color, node_type, status, shape?,
                        entity_type?, parent?, parent_protein? }, … ],
     edges        : [ { source, target, id, interaction, edge_type, category,
                        color, line_style, arrow_shape, width, tooltip,
                        evidence_url? }, … ],
     layout       : { name, rankDir, … },          // dagre options
     container_id : "network-cy",                  // ignored – we use el
     node_font_size : 12
   }
   ========================================================================== */

HTMLWidgets.widget({

  name: "cytoscapeNetwork",
  type: "output",

  /* ── factory ──────────────────────────────────────────────────────────── */
  factory: function (el, width, height) {

    // State kept between renderValue calls so we can destroy cleanly
    var cy      = null;
    var tooltip = null;

    /* helper – open a URL safely in a new tab */
    function openSafe(url) {
      if (!url || typeof url !== "string") return;
      url = url.trim();
      if (!url || url === "NA") return;
      if (!/^https?:\/\//i.test(url)) return;
      var w = window.open(url, "_blank", "noopener,noreferrer");
      if (w) w.opener = null;
    }

    /* helper – build Cytoscape stylesheet from x.node_font_size */
    function buildStyle(nodeFontSize) {
      return [
        /* ── proteins / default nodes ─────────────────────────────────── */
        {
          selector: "node[node_type = 'protein']",
          style: {
            "background-color": "data(color)",
            "label":            "data(label)",
            "shape":            "data(shape)",
            "font-size":        (nodeFontSize || 12) + "px",
            "font-weight":      "bold",
            "color":            "#000",
            "text-valign":      "center",
            "text-halign":      "center",
            "text-wrap":        "wrap",
            "text-max-width":   "140px",
            "border-width":     2,
            "border-color":     "#333",
            "padding":          "5px",
            /* dynamic width/height via mappers */
            "width":  "data(width)",
            "height": "data(height)"
          }
        },
        /* ── PTM child nodes ─────────────────────────────────────────── */
        {
          selector: "node[node_type = 'ptm']",
          style: {
            "shape":            "ellipse",
            "width":            35,
            "height":           35,
            "background-color": "data(color)",
            "border-color":     "#333",
            "border-width":     1.5,
            "label":            "data(label)",
            "font-size":        "10px",
            "font-weight":      "normal",
            "color":            "#000",
            "text-valign":      "center",
            "text-halign":      "center",
            "text-wrap":        "wrap",
            "text-max-width":   "18px"
          }
        },
        /* ── node status (see .node_display_status in R) ───────────────── */
        {
          /* not in the input data: no fill, dashed grey border */
          selector: "node[status = 'latent']",
          style: {
            "background-opacity": 0,
            "border-style":       "dashed",
            "border-color":       "#999",
            "color":              "#555"
          }
        },
        {
          /* in the input, but no logFC to colour by: no fill, solid border */
          selector: "node[status = 'sites_only'], node[status = 'no_logfc']",
          style: {
            "background-opacity": 0,
            "border-style":       "solid"
          }
        },
        {
          /* in the input, but not part of the query: faded */
          selector: "node[status = 'not_queried']",
          style: {
            "background-opacity": 0.4,
            "border-opacity":     0.5
          }
        },
        /* ── invisible compound containers ──────────────────────────── */
        {
          selector: "node[node_type = 'compound']",
          style: {
            "background-opacity": 0,
            "border-width":       0,
            "border-opacity":     0,
            "padding":            "10px",
            "label":              "",
            "z-index":            0
          }
        },
        /* ── all edges (defaults) ────────────────────────────────────── */
        {
          selector: "edge",
          style: {
            "width":                "data(width)",
            "line-color":           "data(color)",
            "line-style":           "data(line_style)",
            "label":                "data(interaction)",
            "curve-style":          "bezier",
            "target-arrow-shape":   "data(arrow_shape)",
            "target-arrow-color":   "data(color)",
            "source-arrow-shape":   "none",
            "source-arrow-color":   "data(color)",
            "edge-text-rotation":   "autorotate",
            "text-margin-y":        -12,
            "text-halign":          "center",
            "font-size":            "11px",
            "font-weight":          "bold",
            "color":                "data(color)",
            "text-background-color":   "#ffffff",
            "text-background-opacity": 0.8,
            "text-background-padding": "2px"
          }
        },
        /* ── undirected (complex) edges ──────────────────────────────── */
        {
          selector: "edge[category = 'complex']",
          style: {
            "line-style":         "solid",
            "target-arrow-shape": "none",
            "source-arrow-shape": "none"
          }
        },
        /* ── phosphorylation edges ───────────────────────────────────── */
        {
          selector: "edge[category = 'phosphorylation']",
          style: {
            "line-style": "dashed",
            "width":      2
          }
        },
        /* ── PTM attachment edges ────────────────────────────────────── */
        {
          selector: "edge[edge_type = 'ptm_attachment']",
          style: {
            "line-style":         "dotted",
            "line-color":         "#9932CC",
            "width":              1.5,
            "target-arrow-shape": "none",
            "source-arrow-shape": "none",
            "label":              "",
            "z-index":            0
          }
        },
        /* ── selected highlight ──────────────────────────────────────── */
        {
          selector: ":selected",
          style: {
            "border-width": 4,
            "border-color": "#FFD700",
            "line-color":   "#FFD700"
          }
        }
      ];
    }

    /* helper – reposition PTM nodes in a small arc below their parent */
    function repositionPTMNodes(cyInstance) {
      var ptmNodes = cyInstance.nodes('[node_type = "ptm"]');
      ptmNodes.forEach(function (ptmNode) {
        var parentId   = ptmNode.data("parent_protein");
        var parentNode = cyInstance.getElementById(parentId);
        if (!parentNode || parentNode.length === 0) return;

        var parentPos = parentNode.position();
        var parentW   = parentNode.outerWidth();
        var parentH   = parentNode.outerHeight();
        var ptmR      = ptmNode.outerWidth() / 2;

        var siblings = cyInstance.nodes('[parent_protein = "' + parentId + '"]');
        var idx      = siblings.indexOf(ptmNode);
        var total    = siblings.length;

        var angleStart = Math.PI * 0.15;
        var angleEnd   = Math.PI * 0.85;
        var angle = (total === 1)
          ? Math.PI / 2
          : angleStart + (angleEnd - angleStart) * (idx / (total - 1));

        var offsetX = (parentW / 2 + ptmR + 4) * Math.cos(angle);
        var offsetY = (parentH / 2 + ptmR + 4) * Math.sin(angle);

        ptmNode.position({
          x: parentPos.x + offsetX,
          y: parentPos.y + offsetY
        });
      });
    }

    /* helper - escape text for the legend's innerHTML */
    function escapeHtml(text) {
      return String(text)
        .replace(/&/g, "&amp;").replace(/</g, "&lt;")
        .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
    }

    /* helper - small SVG of a Cytoscape node shape for the legend */
    function shapeIcon(shape) {
      var body;
      switch (shape) {
        case "hexagon":
          body = '<polygon points="5,2 17,2 22,9 17,16 5,16 0,9"/>'; break;
        case "diamond":
          body = '<polygon points="11,1 21,9 11,17 1,9"/>'; break;
        case "octagon":
          body = '<polygon points="7,1 15,1 21,6 21,12 15,17 7,17 1,12 1,6"/>'; break;
        case "barrel":
          body = '<path d="M4,1 H18 Q22,9 18,17 H4 Q0,9 4,1 Z"/>'; break;
        default:
          body = '<rect x="1" y="2" width="20" height="14" rx="4"/>';
      }
      return '<svg width="22" height="18" style="margin-right:7px;flex-shrink:0;" ' +
             'fill="#fff" stroke="#333" stroke-width="1.5">' + body + '</svg>';
    }

    /* helper - one legend row: a swatch followed by a label */
    function legendRow(swatchHtml, label) {
      return '<div style="display:flex;align-items:center;margin-bottom:5px;font-size:12px;">' +
             swatchHtml + '<span>' + escapeHtml(label) + '</span></div>';
    }

    /* helper - legend heading */
    function legendHeading(text) {
      return '<div style="font-weight:bold;margin:10px 0 6px 0;font-size:13px;">' +
             text + '</div>';
    }

    /* helper – build the legend panel beside the network */
    function buildLegend(cyInstance, legendEl) {
      if (!legendEl) return;

      /* Edge types: one entry per statement type drawn, styled as drawn */
      var edgeTypes = {};
      cyInstance.edges().forEach(function (e) {
        if (e.data("edge_type") === "ptm_attachment") return;
        var type = e.data("interaction") || "";
        if (!type || edgeTypes[type]) return;
        edgeTypes[type] = { color: e.data("color"), style: e.data("line_style") };
      });
      var edgeItems = Object.keys(edgeTypes).sort().map(function (type) {
        var c = edgeTypes[type];
        var line = (c.style === "dashed" || c.style === "dotted")
          ? "border-top:2px " + c.style + " " + c.color + ";"
          : "background-color:" + c.color + ";";
        return legendRow(
          '<div style="width:28px;height:3px;' + line + 'margin-right:7px;flex-shrink:0;"></div>',
          type.replace(/([a-z])([A-Z])/g, "$1 $2"));
      }).join("");

      /* Node status: only the statuses present, worded relative to the input */
      var statusConfigs = [
        { status: "latent",      label: "not in input data",
          box: "border:2px dashed #999;background:transparent;" },
        { status: "sites_only",  label: "no protein-level row in input",
          box: "border:2px solid #333;background:transparent;" },
        { status: "no_logfc",    label: "in input, no logFC shown",
          box: "border:2px solid #333;background:transparent;" },
        { status: "not_queried", label: "in input, not in query",
          box: "border:2px solid rgba(51,51,51,.5);background:rgba(255,165,144,.4);" }
      ];
      var statuses = {};
      cyInstance.nodes().forEach(function (n) {
        if (n.data("status")) statuses[n.data("status")] = true;
      });
      var statusItems = statusConfigs
        .filter(function (c) { return statuses[c.status]; })
        .map(function (c) {
          return legendRow(
            '<div style="width:22px;height:14px;border-radius:4px;box-sizing:border-box;' +
            c.box + 'margin-right:7px;flex-shrink:0;"></div>', c.label);
        })
        .join("");

      /* Node shape: listed when the network has more than one shape */
      var shapeTypes = {};
      cyInstance.nodes('[node_type = "protein"]').forEach(function (n) {
        var shape = n.data("shape") || "round-rectangle";
        var type  = n.data("entity_type") || "";
        if (type === "ptm_site") type = "protein";
        shapeTypes[shape] = shapeTypes[shape] || {};
        if (type) shapeTypes[shape][type] = true;
      });
      var shapes = Object.keys(shapeTypes);
      var shapeItems = shapes.length < 2 ? "" : shapes.map(function (shape) {
        var types = Object.keys(shapeTypes[shape]).sort();
        return legendRow(shapeIcon(shape), types.length ? types.join(", ") : "other");
      }).join("");

      var hasPtm = cyInstance.nodes('[node_type = "ptm"]').length > 0;

      legendEl.innerHTML =
        '<div style="font-weight:bold;margin-bottom:8px;font-size:13px;">Node color (logFC)</div>' +
        '<div style="display:flex;align-items:flex-start;margin-bottom:12px;">' +
        '  <div style="width:18px;height:110px;background:linear-gradient(to top,#ADD8E6,#D3D3D3,#FFA590);border:1px solid #999;border-radius:3px;margin-right:7px;flex-shrink:0;"></div>' +
        '  <div style="display:flex;flex-direction:column;justify-content:space-between;height:110px;font-size:11px;">' +
        '    <span>Upregulated</span><span>Neutral</span><span>Downregulated</span>' +
        '  </div></div>' +
        statusItems +
        (shapeItems ? legendHeading("Node shape") + shapeItems : '') +
        (edgeItems ? legendHeading("Edge types") + edgeItems : '') +
        (hasPtm
          ? '<div style="margin-top:12px;padding:7px;background:#e3f2fd;border-radius:4px;font-size:10px;line-height:1.4;">' +
            '<strong>PTM info:</strong> Hover over edges to see overlapping PTM sites.</div>'
          : '') +
        '<div style="margin-top:8px;padding:7px;background:#fff3cd;border-radius:4px;font-size:10px;line-height:1.4;">' +
        '<strong>Delete edge:</strong> Right-click or Ctrl+Click an edge to remove it from the network.</div>';
    }
    
    // Helper to delete an edge and notify Shiny
    function deleteEdge(edge) {
      if (window.Shiny) {
        Shiny.setInputValue(el.id + "_edge_deleted", {
          source:      edge.data("source"),
          target:      edge.data("target"),
          interaction: edge.data("interaction")
        }, { priority: "event" });
      }
      edge.remove();
    }

    /* ── renderValue ──────────────────────────────────────────────────── */
    return {
      renderValue: function (x) {

        /* Destroy previous instance */
        if (cy) { cy.destroy(); cy = null; }
        if (tooltip) { tooltip.parentNode && tooltip.parentNode.removeChild(tooltip); tooltip = null; }

        el.innerHTML = "";
        var PANEL_W = 160;
        var elH     = el.offsetHeight || height || 600;

        el.style.cssText = [
          "display:flex",
          "width:100%",
          "height:" + elH + "px",
          "box-sizing:border-box"
        ].join(";");

        /* Left: Cytoscape canvas — explicit px so Cytoscape always gets
           real dimensions regardless of flex/CSS resolution order */
        var cyContainer = document.createElement("div");
        cyContainer.style.cssText = [
          "flex:1",
          "min-width:0",
          "height:" + elH + "px"
        ].join(";");

        /* Right panel — shared background for button + legend */
        var PANEL_BG = "#f8f9fa";
        var rightPanel = document.createElement("div");
        rightPanel.style.cssText = [
          "width:" + PANEL_W + "px",
          "flex-shrink:0",
          "display:flex",
          "flex-direction:column",
          "background:" + PANEL_BG,
          "border-left:1px solid #dee2e6",
          "box-sizing:border-box"
        ].join(";");

        /* Button bar — right-aligned inside the panel */
        var btnBar = document.createElement("div");
        btnBar.style.cssText = [
          "display:flex",
          "justify-content:flex-start",
          "padding:8px 10px 6px 10px",
          "background:" + PANEL_BG,
          "border-bottom:1px solid #dee2e6"
        ].join(";");

        var btn = document.createElement("button");
        btn.textContent = "Export PNG";
        btn.style.cssText = [
          "padding:4px 10px",
          "cursor:pointer",
          "font-size:13px",
          "background:#28a745",
          "color:white",
          "border:none",
          "border-radius:4px",
          "font-family:Arial,sans-serif",
          "white-space:nowrap"
        ].join(";");

        btn.addEventListener("click", function () {
          /* Network PNG via Cytoscape */
          var networkPng = cy.png({ output: "base64uri", bg: "white", full: true, scale: 8 });
          var a1 = document.createElement("a");
          a1.href     = networkPng;
          a1.download = "network.png";
          a1.click();

          /* Legend PNG via html2canvas (if available) */
          setTimeout(function () {
            if (typeof html2canvas === "function") {
              html2canvas(legendPanel, { backgroundColor: PANEL_BG, scale: 8 })
                .then(function (canvas) {
                  var a2 = document.createElement("a");
                  a2.href     = canvas.toDataURL("image/png");
                  a2.download = "network_legend.png";
                  a2.click();
                });
            }
          }, 300);
        });

        btnBar.appendChild(btn);

        /* Legend panel — fills remaining vertical space, scrolls if needed */
        var legendPanel = document.createElement("div");
        legendPanel.className = "cytoscape-network-legend";
        legendPanel.style.cssText = [
          "flex:1",
          "overflow-y:auto",
          "padding:10px",
          "font-family:Arial,sans-serif",
          "box-sizing:border-box",
          "background:" + PANEL_BG
        ].join(";");

        rightPanel.appendChild(btnBar);
        rightPanel.appendChild(legendPanel);

        el.appendChild(cyContainer);
        el.appendChild(rightPanel);

        /* Build combined elements array from pre-serialised strings.
           R passes them as an array of JSON-string fragments; we re-parse. */
        var elements = [];
        (x.elements || []).forEach(function (frag) {
          if (typeof frag === "string") {
            try {
              elements.push(JSON.parse(frag));
            } catch (err) {
              console.warn("Skipping invalid element JSON fragment:", err);
            }
          } else if (frag && typeof frag === "object") {
            elements.push(frag);
          }
        });

        /* Layout – merge defaults with whatever R sends */
        var layout = Object.assign({
          name:          "dagre",
          rankDir:       "TB",
          animate:       true,
          fit:           true,
          padding:       30,
          spacingFactor: 1.5,
          nodeSep:       50,
          edgeSep:       20,
          rankSep:       80
        }, x.layout || {});

        /* Initialise Cytoscape */
        cytoscape.use(cytoscapeDagre);   // register dagre layout

        cy = cytoscape({
          container: cyContainer,
          elements:  elements,
          style:     buildStyle(x.node_font_size),
          layout:    layout
        });

        /* After layout, fan PTM nodes around their parent protein */
        cy.on("layoutstop", function () {
          repositionPTMNodes(cy);
          buildLegend(cy, legendPanel);
        });

        /* ── Tooltip ─────────────────────────────────────────────────── */
        tooltip = document.createElement("div");
        tooltip.style.cssText = [
          "position:fixed",
          "background:rgba(0,0,0,.88)",
          "color:#fff",
          "padding:7px 11px",
          "border-radius:4px",
          "font-size:12px",
          "font-family:Arial,sans-serif",
          "pointer-events:none",
          "z-index:99999",
          "box-shadow:0 2px 8px rgba(0,0,0,.3)",
          "display:none",
          "max-width:300px",
          "white-space:pre-wrap",
          "word-wrap:break-word"
        ].join(";");
        document.body.appendChild(tooltip);

        cy.on("mouseover", "edge", function (evt) {
          var txt = evt.target.data("tooltip");
          if (txt && txt.trim()) {
            tooltip.textContent  = txt;
            tooltip.style.display = "block";
          }
        });

        cy.on("mousemove", "edge", function (evt) {
          if (tooltip.style.display === "block") {
            tooltip.style.left = (evt.originalEvent.clientX + 12) + "px";
            tooltip.style.top  = (evt.originalEvent.clientY - 28) + "px";
          }
        });

        cy.on("mouseout", "edge", function () {
          tooltip.style.display = "none";
        });

        /* ── Edge tap: Ctrl+Click → delete; plain click → evidence link ── */
        cy.on("cxttap tap", "edge", function (evt) {
          var edge = evt.target;
          // skip ptm attachment edges
          if (edge.data("edge_type") === "ptm_attachment") return;

          // Ctrl+Click or Right Click → delete edge
          if (evt.type === "cxttap" || (evt.originalEvent && evt.originalEvent.ctrlKey)) {
            deleteEdge(edge);
            buildLegend(cy, legendPanel);
            return;
          }

          // Plain click → open evidence link
          openSafe(edge.data("evidence_url"));
          if (window.Shiny) {
            Shiny.setInputValue(el.id + "_edge_clicked", {
              source:       edge.data("source"),
              target:       edge.data("target"),
              interaction:  edge.data("interaction"),
              edge_type:    edge.data("edge_type"),
              category:     edge.data("category"),
              evidence_url: edge.data("evidence_url")
            });
          }
        });

        /* ── Node click — report to Shiny ───────────────────────────── */
        cy.on("tap", "node", function (evt) {
          var node = evt.target;
          // skip compound and ptm satellite nodes
          if (node.data("node_type") === "compound" || 
            node.data("node_type") === "ptm") return;
          if (window.Shiny) {
            Shiny.setInputValue(el.id + "_node_clicked", {
              id:        node.data("id"),
              label:     node.data("label"),
              color:     node.data("color"),
              node_type: node.data("node_type")
            });
          }
        });

        /* ── Expose cy instance for external access (e.g. Shiny) ─────── */
        el._cytoscapeInstance = cy;

        /* ── Build legend in sibling element (if present) ────────────── */
        var legendEl = el.parentNode
          ? el.parentNode.querySelector(".cytoscape-network-legend")
          : null;
        buildLegend(cy, legendEl);
      },

      /* ── resize ───────────────────────────────────────────────────── */
      resize: function (newWidth, newHeight) {
        if (cy) {
          cy.resize();
          cy.fit();
        }
      }
    };
  }
});
