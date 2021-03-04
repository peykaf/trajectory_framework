import json
import pickle
from scipy.interpolate import interp1d

class beamWeight:

    def __init__(self, path):
        self.path = path
        self.get_beam_weight()

    def get_beam_weight(self):
        with open(self.path) as f:
            data = json.load(f)
        
        self.beam_weights = []
        for cpt in data['control_points']:
            self.beam_weights.append(sum(cpt['beamlet_weights']))

        min_val = min(self.beam_weights)
        max_val = max(self.beam_weights)
        m = interp1d([min_val, max_val],[0.05, 0.95])
        self.beam_weights = [float(m(self.beam_weights[i])) for i in range(len(self.beam_weights))]
        
        with open('score/fmo.pkl', 'wb') as f:
            pickle.dump(self.beam_weights, f)

if __name__ == "__main__":
    data_path = 'build/combined_doses/dec11_4.weights'
    fmo = beamWeight(data_path)
            