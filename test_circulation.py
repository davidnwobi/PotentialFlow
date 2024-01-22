from src import circulation as c
from src import elementary_flows as ef
from src import data_collections as dc


def test_circulation_uniform_flow():
    '''
    Tests the circulation of an ellipse placed in a uniform flow field.
    '''
    # Define the flow field
    u1 = ef.UniformFlow(horizontal_vel=1, vertical_vel=0)

    # Define the ellipse
    ellipse_def = dc.Ellipse(x0=0, y0=0, a=1, b=1)

    # Compute the ellipse properties
    ellipse_properties = c.compute_ellipse_and_circulation(u1, ellipse_def)

    # Assert the circulation is zero
    assert ellipse_properties.circulation == 0

def test_circulation_non_lifting_flow():
    '''
    Tests the circulation of an ellipse placed in a non-lifting flow field.
    '''
    # Define the flow field
    u1 = ef.UniformFlow(horizontal_vel=1, vertical_vel=0)
    d1 = ef.Doublet(x=0, y=0, kappa=1)
