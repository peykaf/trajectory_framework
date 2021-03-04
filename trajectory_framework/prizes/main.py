from utils.options import get_options
from beams import get_beams

def run(opts):
    phi, theta = get_beams.full_path(opts)

if __name__ == "__main__":
    run(get_options())