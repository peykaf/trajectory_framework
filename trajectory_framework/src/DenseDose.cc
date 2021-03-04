#include <vector>
#include <string>
#include "json.hh"
#include "DenseDose.hh"

void DenseDose::from_3ddose(std::string filename, bool norm)
{
    std::string buffer;
    char *token;
    // Save the input filename to know how to name the final dose.
    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer); // number of voxels
    std::getline(beamlet_file, buffer); // x voxel coordinates
    std::getline(beamlet_file, buffer); // y voxel coordinates
    std::getline(beamlet_file, buffer); // z voxel coordinates
    std::getline(beamlet_file, buffer); // Dose values

    char *buf_dup = strdup(buffer.c_str());
    double max_dose = 0.0;

    token = strtok(buf_dup, " ");
    int idx = 0;
    while (token != NULL)
    {
        grid.push_back(atof(token));
        if (grid[idx] > max_dose)
        {
            max_dose = grid[idx];
        }
        idx++;
        token = strtok(NULL, " ");
    }

    if (norm)
    {
        // Ideally I really don't want to be doing this...
        double norm = 1.0 / max_dose;

        for (size_t i = 0; i < grid.size(); i++)
        {
            grid[i] = grid[i] * norm;
        }
    }

    free(buf_dup);
    beamlet_file.close();
}

void DenseDose::from_3ddose(std::string filename, size_t total_vox, bool norm)
{
    std::string buffer;
    // Save the input filename to know how to name the final dose.
    std::ifstream beamlet_file(filename);
    std::getline(beamlet_file, buffer); // number of voxels
    std::getline(beamlet_file, buffer); // x voxel coordinates
    std::getline(beamlet_file, buffer); // y voxel coordinates
    std::getline(beamlet_file, buffer); // z voxel coordinates
    std::getline(beamlet_file, buffer); // Dose values

    grid.resize(total_vox);

    //char *buf_dup = strdup(buffer.c_str());
    char *buf_dup = new char[buffer.length() + 1];
    memcpy(buf_dup, buffer.c_str(), buffer.length() + 1);
    double max_dose = 0.0;

    // This is really not pretty... but it's faster than stringstreams
    // and Boost's tokenizer.
    grid[0] = atof(strtok(buf_dup, " "));
    for (size_t i = 1; i < total_vox; i++)
    {
        grid[i] = atof(strtok(NULL, " "));
        if (grid[i] > max_dose)
        {
            max_dose = grid[i];
        }
    }

    if (norm)
    {
        // Ideally I really don't want to be doing this...
        double norm = 1.0 / max_dose;

        for (size_t i = 0; i < total_vox; i++)
        {
            grid[i] = grid[i] * norm;
        }
    }

    //free(buf_dup);
    beamlet_file.close();
}

void DenseDose::from_minidos(std::string filename, bool norm)
{
    int num_voxels[3], num_nonzero;
    float voxel_size[3];
    float topleft[3];
    size_t total_vox;

    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);
    beamlet_file.read((char *)num_voxels, 3 * sizeof(int));
    beamlet_file.read((char *)voxel_size, 3 * sizeof(float));
    beamlet_file.read((char *)topleft, 3 * sizeof(float));
    beamlet_file.read((char *)&num_nonzero, sizeof(int));

    total_vox = num_voxels[0] * num_voxels[1] * num_voxels[2];
    grid = std::vector<double>(total_vox, 0.0);
    //grid.resize(total_vox);
    //std::fill(grid.begin(), grid.end(), 0.0);
    float max_dose = 0.0;
    int vox;
    float dose_val;
    for (int i = 0; i < num_nonzero; i++)
    {
        beamlet_file.read((char *)&vox, sizeof(int));
        beamlet_file.read((char *)&dose_val, sizeof(float));
        grid[vox] = dose_val;
        if (dose_val > max_dose)
            max_dose = dose_val;
    }

    if (norm)
    {
        // Ideally I really don't want to be doing this...
        double norm = 1.0 / max_dose;

        for (size_t i = 0; i < grid.size(); i++)
        {
            grid[i] = grid[i] * norm;
        }
    }

    beamlet_file.close();
}

void DenseDose::from_bindos(std::string filename, bool norm)
{
    int num_voxels[3], num_nonzero;
    size_t total_vox;


    std::ifstream beamlet_file(filename, std::ios::in | std::ios::binary);

    if (!beamlet_file.is_open())
    {
        std::cout << "Could not find file: " << filename << std::endl;
        throw;
    }
    beamlet_file.read((char *)num_voxels, 3 * sizeof(int));

    float *x_voxels = new float[num_voxels[0] + 1];
    float *y_voxels = new float[num_voxels[1] + 1];
    float *z_voxels = new float[num_voxels[2] + 1];

    beamlet_file.read((char *)x_voxels, (num_voxels[0] + 1) * sizeof(float));
    beamlet_file.read((char *)y_voxels, (num_voxels[1] + 1) * sizeof(float));
    beamlet_file.read((char *)z_voxels, (num_voxels[2] + 1) * sizeof(float));

    beamlet_file.read((char *)&num_nonzero, sizeof(int));

    int *voxels = new int[num_nonzero];
    float *doses = new float[num_nonzero];

    beamlet_file.read((char *)voxels, num_nonzero * sizeof(int));
    beamlet_file.read((char *)doses, num_nonzero * sizeof(float));

    float max_dose = 0.0;
    total_vox = num_voxels[0] * num_voxels[1] * num_voxels[2];
    grid = std::vector<double>(total_vox, 0.0);
    for (size_t i = 0; i < static_cast<size_t>(num_nonzero); i++)
    {
        grid[voxels[i]] = doses[i];
        if (doses[i] > max_dose)
            max_dose = doses[i];
    }

    beamlet_file.close();

    if (norm)
    {
        // Ideally I really don't want to be doing this...
        double norm = 1.0 / max_dose;

        for (size_t i = 0; i < grid.size(); i++)
        {
            grid[i] = grid[i] * norm;
        }
    }

    delete[] x_voxels;
    delete[] y_voxels;
    delete[] z_voxels;
    delete[] voxels;
    delete[] doses;
}

void DenseDose::from_file(std::string filename, bool norm)
{
    std::string extension = find_extension(filename);
    if (extension == "3ddose")
    {
        this->from_3ddose(filename, norm);
    }
    else if (extension == "minidos")
    {
        this->from_minidos(filename, norm);
    }
    else if (extension == "bindos")
    {
        this->from_bindos(filename, norm);
    }
    else
    {
        throw "Dose file extension not recognised";
    }
}

void DenseDose::from_file(std::string filename, size_t total_vox, bool norm)
{
    std::string extension = find_extension(filename);
    if (extension == "3ddose")
    {
        this->from_3ddose(filename, total_vox, norm);
    }
    else if (extension == "minidos")
    {
        this->from_minidos(filename, norm);
    }
    else if (extension == "bindos")
    {
        this->from_bindos(filename, norm);
    }
    else
    {
        throw "Dose file extension not recognised";
    }
}
