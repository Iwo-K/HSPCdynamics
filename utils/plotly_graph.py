import networkx as nx
import plotly.graph_objects as go
import numpy as np
import math


def draw_graph_pl(G,
                  node_positions='pos',
                  # node_names='leiden',
                  node_colors='ln(self-renewal<sup>-1</sup>)',
                  node_sizes='leiden_sizes',
                  node_hovertext='text',
                  edge_widths='flux',
                  edge_width_scale=0.36,
                  edge_hovertext='text'):

    node_x = []
    node_y = []
    colors = []
    sizes = []
    for node in G.nodes():
        x, y = G.nodes[node][node_positions]
        node_x.append(x)
        node_y.append(y)
        col = G.nodes[node][node_colors]
        colors.append(col)
        s = (G.nodes[node][node_sizes]**0.25)*25
        sizes.append(s)

    # Creating nodes (with traces)
    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale='YlGnBu',
            # reversescale=True,
            colorscale='plasma',
            color=colors,
            size=sizes,
            colorbar=dict(
                thickness=15,
                title=node_colors,
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_trace.text = list(nx.get_node_attributes(G, node_hovertext).values())
    node_trace.marker.color = list(nx.get_node_attributes(G, node_colors).values())


    # Creating edges (added as paths - bezier curves and triangles
    paths = []
    mid_x = []
    mid_y = []
    for edge in G.edges():
        a0 = np.array(G.nodes[edge[0]][node_positions])
        a1 = np.array(G.nodes[edge[1]][node_positions])
        a0, a1 = shrink_points(a0, a1, step=0.90)

        vec = np.array(a1-a0)
        l = np.linalg.norm(vec)
        vec90 = np.array([vec[1], -vec[0]])
        ac = np.array([
            a0[0] + 0.5*vec[0] + vec90[0]/5,
            a0[1] + 0.5*vec[1] + vec90[1]/5
        ])


        p = f"M {a0[0]},{a0[1]} Q {ac[0]},{ac[1]} {a1[0]},{a1[1]}"

        w = G.edges()[edge][edge_widths]
        p = dict(type='path',
                 path=p,
                 line_color='Grey',
                 line=dict(width=w*edge_width_scale))
        if w < 0.05:
            p['line']['width'] = 0.4
            p['line']['dash'] = 'dash'
        paths.append(p)

        # Adding arrows
        arrow = make_arrowhead(tip=a1, vec=vec, size=np.log(5+w)*0.1)
        paths.append(arrow)

        # Storing coordinates for invisible edge midpoints
        mid = np.array([
            a0[0] + 0.5*vec[0] + vec90[0]/10,
            a0[1] + 0.5*vec[1] + vec90[1]/10
        ])
        mid_x.append(mid[0])
        mid_y.append(mid[1])

    # Creating invisible edge midpoints (with traces) for the edges to permit hovertext

    mid_trace = go.Scatter(
        x=mid_x,
        y=mid_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            opacity=0,
            color='Green',
            size=20))

    mid_trace.text = list(nx.get_edge_attributes(G, edge_hovertext).values())

    # Assembling the figure
    fig = go.Figure(data=[node_trace, mid_trace],
                    layout=go.Layout(
                        title='<br>Network graph made with Python',
                        titlefont_size=16,
                        plot_bgcolor='white',
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        xaxis=dict(showgrid=True, zeroline=False, showticklabels=False, fixedrange=True),
                        yaxis=dict(showgrid=True, zeroline=False, showticklabels=False, fixedrange=True))
                )
    fig.update_layout(shapes=paths)

    fig.show()
    return(fig)


def make_arrowhead(tip, vec, size=0.75, theta=-0.28, arrownudge=0.9):
    # Consider using just two lines to create arrowhead, alignment will be easier
    # and will scale appriopriately
    l = np.linalg.norm(vec)
    vec = vec*size/l

    #Rotating
    c, s = np.cos(theta), np.sin(theta)
    R = np.array(((c, -s), (s, c)))
    vec = np.matmul(vec, R)

    vec90 = np.array([vec[1], -vec[0]])

    tip = tip + arrownudge*vec
    v1 = tip
    v2 = tip - vec + vec90
    v3 = tip - vec - vec90
    out = dict(type='path',
               path=f"M {v1[0]} {v1[1]} L {v2[0]} {v2[1]} L {v3[0]} {v3[1]} Z",
               fillcolor='Grey',
               line_color='Grey')
    return(out)

def shrink_points(x, y, step = 0.5):
    vec = y-x
    l = np.linalg.norm(vec)

    x = x+step*vec/l
    y = y-step*vec/l
    return(x, y)


    # annots = []
    #     print(a1)
    #     print(vec)
    #     annot = dict(ax=a1[0]- vec[0], ay=a1[1] - vec[1],
    #                  x=a1[0], y=a1[1],
    #                  axref='x', ayref='y',
    #                  xref='x', yref='y',
    #                  showarrow=True,
    #                  arrowhead=1,
    #                  arrowsize=50,
    #                  arrowwidth=0.1)
    #     annots.append(annot)
