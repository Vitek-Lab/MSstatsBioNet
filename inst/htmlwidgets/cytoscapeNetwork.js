/* ==========================================================================
   cytoscapeNetwork.js
   htmlwidgets binding for the cytoscapeNetwork R package.

   Dependencies (declared in cytoscapeNetwork.yaml):
     - cytoscape.min.js
     - dagre.min.js / graphlib.min.js
     - cytoscape-dagre.js

   The `x` object passed from R (via createWidget / jsonlite serialisation):
   {
     nodes        : [ { id, label, color, node_type, parent?, parent_protein? }, … ],
     edges        : [ { source, target, id, interaction, edge_type, category,
                        color, line_style, arrow_shape, width, tooltip,
                        evidenceLink? }, … ],
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
            "shape":            "round-rectangle",
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
            "width":  "mapData(label.length, 0, 20, 60, 150)",
            "height": 40
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
        /* ── bidirectional edges – source arrow too ──────────────────── */
        {
          selector: "edge[edge_type = 'bidirectional']",
          style: {
            "source-arrow-shape": "triangle",
            "target-arrow-shape": "triangle"
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

    /* helper – build the legend panel beside the network */
    function buildLegend(cyInstance, legendEl) {
      if (!legendEl) return;

      var edgeTypeConfigs = [
        { type: "Activation",     color: "#44AA44", label: "Activation",      dash: false },
        { type: "Inhibition",     color: "#FF4444", label: "Inhibition",      dash: false },
        { type: "IncreaseAmount", color: "#4488FF", label: "Increase Amount", dash: false },
        { type: "DecreaseAmount", color: "#FF8844", label: "Decrease Amount", dash: false },
        { type: "Phosphorylation",color: "#9932CC", label: "Phosphorylation", dash: true  },
        { type: "Complex",        color: "#8B4513", label: "Complex",         dash: false }
      ];

      var existingTypes = {};
      cyInstance.edges().forEach(function (e) {
        var raw = e.data("interaction") || "";
        existingTypes[raw.replace(" (bidirectional)", "")] = true;
      });

      var edgeItems = edgeTypeConfigs
        .filter(function (c) { return existingTypes[c.type]; })
        .map(function (c) {
          var dash = c.dash ? "border-top: 2px dashed " + c.color + ";" : "background-color:" + c.color + ";";
          return '<div style="display:flex;align-items:center;margin-bottom:5px;font-size:12px;">' +
                 '<div style="width:28px;height:3px;' + dash + 'margin-right:7px;flex-shrink:0;"></div>' +
                 '<span>' + c.label + '</span></div>';
        })
        .join("");

      legendEl.innerHTML =
        '<div style="font-weight:bold;margin-bottom:8px;font-size:13px;">Node color (logFC)</div>' +
        '<div style="display:flex;align-items:flex-start;margin-bottom:12px;">' +
        '  <div style="width:18px;height:110px;background:linear-gradient(to top,#ADD8E6,#D3D3D3,#FFA590);border:1px solid #999;border-radius:3px;margin-right:7px;flex-shrink:0;"></div>' +
        '  <div style="display:flex;flex-direction:column;justify-content:space-between;height:110px;font-size:11px;">' +
        '    <span>Upregulated</span><span>Neutral</span><span>Downregulated</span>' +
        '  </div></div>' +
        (edgeItems ? '<div style="font-weight:bold;margin-bottom:6px;font-size:13px;">Edge types</div>' + edgeItems : '') +
        '<div style="margin-top:12px;padding:7px;background:#e3f2fd;border-radius:4px;font-size:10px;line-height:1.4;">' +
        '<strong>PTM info:</strong> Hover over edges to see overlapping PTM sites.</div>';
    }

    /* ── renderValue ──────────────────────────────────────────────────── */
    return {
      renderValue: function (x) {

        /* Destroy previous instance */
        if (cy) { cy.destroy(); cy = null; }
        if (tooltip) { tooltip.parentNode && tooltip.parentNode.removeChild(tooltip); tooltip = null; }

        /* Ensure the container has explicit pixel dimensions */
        el.style.width  = el.style.width  || width  + "px";
        el.style.height = el.style.height || height + "px";

        /* Build combined elements array from pre-serialised strings.
           R passes them as an array of JSON-string fragments; we re-parse. */
        var elements = (x.elements || []).map(function (frag) {
          return (typeof frag === "string") ? JSON.parse(frag) : frag;
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
        
        el.innerHTML = "";  // clear on re-render
        // Outer flex wrapper — fills the widget element
        var wrapper = document.createElement("div");
        wrapper.style.cssText = "display:flex;width:100%;height:100%;";
        
        // Left: Cytoscape canvas
        var cyContainer = document.createElement("div");
        cyContainer.style.cssText = "flex:1;height:100%;min-width:0;";
        
        // Right: legend panel
        var legendPanel = document.createElement("div");
        legendPanel.className = "cytoscape-network-legend";
        legendPanel.style.cssText = [
          "width:180px",
          "flex-shrink:0",
          "padding:12px",
          "background:#f8f9fa",
          "border-left:1px solid #dee2e6",
          "overflow-y:auto",
          "font-family:Arial,sans-serif",
          "box-sizing:border-box"
        ].join(";");
        
        wrapper.appendChild(cyContainer);
        wrapper.appendChild(legendPanel);
        el.appendChild(wrapper);

        cy = cytoscape({
          container: cyContainer,
          elements:  elements,
          style:     buildStyle(x.node_font_size),
          layout:    layout
        });
        
        // Inject an export PNG button above the container
        var btnBar = document.createElement("div");
        btnBar.style.cssText = "display:flex;justify-content:flex-end;margin-bottom:6px;";
        
        var btn = document.createElement("button");
        btn.textContent = "Export PNG";
        btn.style.cssText = [
          "padding:5px 12px",
          "cursor:pointer",
          "font-size:12px",
          "background:#28a745",
          "color:white",
          "border:none",
          "border-radius:4px",
          "font-family:Arial,sans-serif"
        ].join(";");
        
        btn.addEventListener("click", function () {
          var png = cy.png({
            output: "base64uri",
            bg:     "white",
            full:   true,
            scale:  3
          });
          var a = document.createElement("a");
          a.href     = png;
          a.download = "network.png";
          a.click();
        });
        
        btnBar.appendChild(btn);
        el.parentNode.insertBefore(btnBar, el);

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

        /* ── Evidence link on edge click ─────────────────────────────── */
        cy.on("tap", "edge", function (evt) {
          var link = evt.target.data("evidenceLink");
          openSafe(link);
        });

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
