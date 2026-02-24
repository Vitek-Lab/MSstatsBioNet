# Generate Cytoscape visualization configuration

This function creates a complete Cytoscape configuration object that can
be used to render a network visualization. It's decoupled from any
specific UI framework.

## Usage

``` r
generateCytoscapeConfig(
  nodes,
  edges,
  display_label_type = "id",
  container_id = "network-cy",
  event_handlers = NULL,
  layout_options = NULL,
  node_font_size = 12
)
```

## Arguments

- nodes:

  List of nodes from getSubnetworkFromIndra

- edges:

  List of edges from getSubnetworkFromIndra

- display_label_type:

  column of nodes table for displaying node names

- container_id:

  ID of the HTML container element (default: 'network-cy')

- event_handlers:

  Optional list of event handler configurations

- layout_options:

  Optional list of layout configuration options

## Value

List containing: - elements: Combined node and edge elements - style:
Cytoscape style configuration - layout: Layout configuration -
container_id: Container element ID - js_code: Complete JavaScript code
(for backward compatibility)
