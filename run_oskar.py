#!/usr/bin/env python3

import oskar
import pprint
import datetime
import argparse

def make_vis_data(base_settings, out_ms_file, sky_model_filename, start_freq, num_chan, num_time_steps, int_time, uv_min, uv_max):
    """Run simulation using supplied settings."""

    settings = oskar.SettingsTree('oskar_sim_interferometer', settings_file=base_settings)

    settings['sky/oskar_sky_model/file'] = sky_model_filename
    settings['observation/start_frequency_hz'] = start_freq * 1e6
    settings['observation/num_channels'] = num_chan
    settings['observation/length'] = str(datetime.timedelta(seconds=num_time_steps * int_time))
    settings['observation/num_time_steps'] = num_time_steps
    settings['interferometer/ms_filename'] = out_ms_file
    settings['interferometer/uv_filter_min'] = uv_min
    settings['interferometer/uv_filter_max'] = uv_max
    settings['interferometer/uv_filter_units'] = 'Metres'
    settings['simulator/use_gpus']=True
    settings['simulator/max_sources_per_chunk'] = 16384
    print('Simulation settings:')
    print(settings.to_dict())

    sim = oskar.Interferometer(settings=settings)
    sim.run()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base_settings')
    parser.add_argument('out_ms_file')
    parser.add_argument('sky_model_filename')
    parser.add_argument('--start_freq', type=float, default=106.0, help='Starting frequency')
    parser.add_argument('--num_chan', type=int, default=1, help='Number of channels')
    parser.add_argument('--num_time_steps', type=int, default=1440, help='Number of time steps')
    parser.add_argument('--int_time', type=int, default=10, help='Integration time')
    parser.add_argument('--uv_min', type=int, default=0, help='Min uv')
    parser.add_argument('--uv_max', type=int, default=1e9, help='Max uv')
    args = parser.parse_args()

    make_vis_data(args.base_settings, args.out_ms_file, args.sky_model_filename,
                  args.start_freq, args.num_chan, args.num_time_steps, args.int_time, args.uv_min, args.uv_max)
