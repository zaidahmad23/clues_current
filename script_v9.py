"""
Automate running of different scripts

python3 automation_script.py -p0 0.5 -s 0.05 -n 5000 --nder 4 --ton 10000 --toff 10500  \
    --converted-filename relate_input \
    --path-to-relate-bin ./relate_v1.1.2_x86_64_static/bin/Relate \
    --relate-output-filename relate_step_1 \
    --relate-map-file-path clues/example/genetic_map.txt \
    --sample-branch-length-coal-file-path clues/example/example.coal \
    --path-to-sample-branch-length-script ./relate_v1.1.2_x86_64_static/scripts/SampleBranchLengths/SampleBranchLengths.sh \
    --sample-branch-length-script-format b \
    --sample-branch-length-script-output-filename clues_input_from_relate \
    --sample-branch-length-script-n-samples 5 \
    --inference-script-output-filename clues_output \
    --output-directory output --create-ancient-samples --case2-script-ancient-samples-generation-gap 500 \
    --case2-script-number-of-ancient-samples 100
"""
import subprocess
import argparse
import os
import sys
import shutil

from pathlib import Path

EXTERNAL_DEPENDENCIES = ["git", "gcc", "Rscript"]


def is_debug_mode_active():
    return os.environ.get("DEBUG", False)

def print_if_debug_mode_active(obj):
    if is_debug_mode_active():
        print(obj)


def execute_command(args, cwd=None):
    try:
        print_if_debug_mode_active(args)
        subprocess.run(args, cwd=cwd, check=True, capture_output=True)
    except subprocess.CalledProcessError as err:
        stderr = err.stderr.decode("utf-8")
        stdout = err.stdout.decode("utf-8")
        if stderr:
            print(stderr)
        elif stdout is not None:
            print(stdout)
        else:
            print(f"Error occurred while executing {args}")
        sys.exit(1)


def clone_rhps_coalescent_repo():
    if not os.path.isdir("rhps_coalescent"):
        execute_command(["git", "clone", "https://github.com/mdedge/rhps_coalescent.git"])


def clone_clues_repo():
    if not os.path.isdir("clues"):
        execute_command(["git", "clone", "https://github.com/35ajstern/clues.git"])


def compile_mssel():
    execute_command(
        ["gcc", "-O2", "-o", "mssel", "mssel.c", "rand1.c", "streecsel.c", "-lm"], cwd="./rhps_coalescent/msseldir"
    )


def run_step(p_initial, s, n, output_file_path, ton, toff):
    command = list(
        map(
            str,
            ['python3', 'step.py', '-p0', p_initial, '-s', s, '-n', n, '--output-file-path', output_file_path, '--ton', ton, '--toff', toff]
        )
    )
    execute_command(command)

def run_case2(p_initial, s, n, ton, toff, ancient_sample_generation_gap, number_of_ancient_samples, output_file_path):
    command = list(
        map(
            str,
            [
                'python3', 'case2.py', '--initial-allele-freq', p_initial, '--selection-coefficient', s,
                '--ton', ton, '--toff', toff, '--effective-population-size', n,
                '--ancient-samples-generation-gap', ancient_sample_generation_gap,
                '--number-of-ancient-samples', number_of_ancient_samples,
                '--output-file-path', output_file_path
            ]
        )
    )
    execute_command(command)


def run_mssel(
    nchroms,
    nreps,
    nder,
    nanc,
    path_to_trajfile,
    sel_spot,
    recombination_rate,
    genome_length,
    mutation_rate,
    output_file,
):
    command = list(
        map(
            str,
            [
                "./rhps_coalescent/msseldir/mssel",
                nchroms,
                nreps,
                nanc,
                nder,
                path_to_trajfile,
                sel_spot,
                "-r",
                recombination_rate,
                genome_length,
                "-t",
                400,
            ],
        )
    )
    with open(output_file, "w") as f:
        try:
            print_if_debug_mode_active(command)
            subprocess.run(command, stdout=f, check=True)
        except subprocess.CalledProcessError as err:
            print(err.stderr if err.stderr else err.stdout)
            sys.exit(1)


def convert_txt_to_haps_and_sample(r_script_path, input_txt_file_path, output_file_name):
    execute_command(f"Rscript {r_script_path} {input_txt_file_path} {output_file_name} 1000000 400".split(" "))


def run_relate(
    relate_binary_path,
    mode,
    haps_file_path,
    sample_file_path,
    map_file_path,
    effective_population_size,
    mutation_rate,
    output_file_path,
):
    command = list(
        map(
            str,
            [
                relate_binary_path,
                "--mode", mode,
                "--haps", haps_file_path,
                "--sample", sample_file_path,
                "--map", map_file_path,
                "-N", effective_population_size,
                "-m", mutation_rate,
                "-o", output_file_path,
            ],
        )
    )
    execute_command(command)


def run_sample_branch_length(
    script_path, input_file_name, mutation_rate, coal_file_path, _format, output_file_name, n_samples
):
    command = list(
        map(
            str,
            [
                script_path,
                "-i", input_file_name,
                "-o", output_file_name,
                "-m", mutation_rate,
                "--coal", coal_file_path,
                "--format", _format,
                "--num_samples", n_samples,
                # "--first_bp", 2500,
                # "--last_bp", 2500
            ],
        )
    )
    execute_command(command)


def run_inference(coalescence_times, ancient_samples_file_path, inference_output_filename):
    command = ["python3", "inference.py"]
    if coalescence_times is not None:
        command += ["--times", coalescence_times]
    if ancient_samples_file_path is not None:
        command += ["--ancientSamps", ancient_samples_file_path]

    command += ["--out", inference_output_filename]
    execute_command(command, cwd="./clues")


def plot(input_file_path, output_file_path):
    execute_command(
        [
            "python3",
            "plot_traj.py",
            input_file_path,
            output_file_path,
        ],
        cwd="./clues",
    )


def fill_defaults(args):
    if args.nanc is None:
        args.nanc = args.nchroms - args.nder


def validate_input_bounds(args):
    if not (0 <= args.nder <= args.nchroms):
        raise argparse.ArgumentTypeError(
            "Argument --nder=%s must greater than or equal to 0 and less than or equal to argument --nchroms=%s"
            % (args.nder, args.nchroms)
        )


    if (
        args.inference_script_coalescence_times_file_path is None and
        args.create_ancient_samples is False
    ):
        raise argparse.ArgumentTypeError(
            "You must specify arguments to at least one of --inference-script-coalescence-times-file-path "
            "or --create-ancient-samples"
        )

    if args.create_ancient_samples:
        if args.case2_script_ancient_samples_generation_gap is None:
            raise argparse.ArgumentTypeError(
                "When specifying --create_ancient_samples, parameters "
                "to --case2-script-ancient-samples-generation-gap are required"
            )
        if args.case2_script_number_of_ancient_samples is None:
            raise argparse.ArgumentTypeError(
                "When specifying --create_ancient_samples, parameters "
                "to --case2-script-number-of-ancient-samples are required"
            )


def ensure_external_dependencies():
    for program in EXTERNAL_DEPENDENCIES:
        if shutil.which(program) is None:
            print(f"{program} is a required dependency. Please install it before running the script.")
            sys.exit(1)

    clone_rhps_coalescent_repo()
    clone_clues_repo()
    compile_mssel()


def ensure_internal_dependencies(args):
    required_files = [
        args.step_script_path,
        args.path_to_converter_script,
        args.path_to_relate_bin,
        args.path_to_sample_branch_length_script,
    ]

    for _file in required_files:
        if not os.path.isfile(_file):
            print("Error: File {args.step_script_path} not found!")


def main():
    # Argument parsing/validation

    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    argparser.add_argument("--output-directory", type=str, default='./')
    argparser.add_argument("--step-script-path", type=str, default="step.py", help="Path to step.py script")
    argparser.add_argument("--step-script-output-file-path", type=str, default="mssel.traj")
    argparser.add_argument(
        "-p0", "--initial-allele-freq", type=float, required=True, help="Present-day allele frequency."
    )
    argparser.add_argument("-s", "--selection-coefficient", type=float, required=True, help="Selection coefficient.")
    argparser.add_argument(
        "-n", "--effective-population-size", type=float, required=True, help="Effective population size."
    )
    argparser.add_argument("--ton", type=int, required=True, help="Time before present that selection starts.")
    argparser.add_argument("--toff", type=int, required=True, help="Time before present that selection ends.")

    argparser.add_argument("--nchroms", type=int, default=400, help="Sample size (in # of chromosomes).")
    argparser.add_argument("--nreps", type=int, default=1, help="Number of independent, identical replicates to run.")
    argparser.add_argument(
        "--nder", type=int, required=True, help="The number of chromosomes that carry the derived allele."
    )
    argparser.add_argument("--nanc", type=int, help="The number of chromosomes that carry the ancestral allele.")

    argparser.add_argument(
        "--genome-length", type=int, default=1000000, help="Length of the genome sequence, in base pairs."
    )
    argparser.add_argument(
        "--sel-spot",
        type=int,
        default=500000,
        help="Position of the allele with the trajectory specified by the trajfile.",
    )
    argparser.add_argument("--mutation-rate", type=float, default=400, help="Population-scaled mutation rate.")
    argparser.add_argument("--recombination-rate", type=int, default=400, help="Population-scaled recombination rate.")
    argparser.add_argument("--massel-output", type=str, default="mssel_output.txt", help="Mssel output filename/path.")

    argparser.add_argument(
        "--path-to-converter-script",
        type=str,
        default="ms2haps_mod.R",
        help="Path/filename of R script that converts .txt to .haps and .sample",
    )
    argparser.add_argument("--converted-filename", type=str, help="Filename of converted .haps and .sample")

    argparser.add_argument("--path-to-relate-bin", type=str, required=True, help="Path to relate binary")
    argparser.add_argument("--relate-mode", type=str, default="All", help="Relate mode")
    argparser.add_argument("--relate-map-file-path", type=str, required=True, help="Relate input map file")
    argparser.add_argument("--relate-output-filename", type=str, required=True, help="Relate output files")


    argparser.add_argument("--path-to-sample-branch-length-script", type=str, required=True)
    argparser.add_argument("--sample-branch-length-script-format", type=str, required=True)
    argparser.add_argument("--sample-branch-length-script-output-filename", type=str, required=True)
    argparser.add_argument("--sample-branch-length-script-n-samples", type=int, required=True)
    argparser.add_argument("--sample-branch-length-coal-file-path", type=str, required=True)

    argparser.add_argument("--inference-script-output-filename", type=str, required=True)
    argparser.add_argument("--inference-script-coalescence-times-file-path", type=str)
    argparser.add_argument("--create-ancient-samples", action="store_true")
    argparser.add_argument("--case2-script-ancient-samples-generation-gap", type=str)
    argparser.add_argument("--case2-script-number-of-ancient-samples", type=int)

    args = argparser.parse_args()
    fill_defaults(args)
    validate_input_bounds(args)

    # Ensuring internal/external dependencies
    ensure_internal_dependencies(args)
    ensure_external_dependencies()

    # Create output directory
    args.output_directory = os.path.abspath(args.output_directory)
    os.makedirs(args.output_directory, exist_ok=True)

    # Update output files
    args.step_script_output_file_path = os.path.join(args.output_directory, args.step_script_output_file_path)
    args.massel_output = os.path.join(args.output_directory, args.massel_output)
    args.converted_filename = os.path.join(args.output_directory, args.converted_filename)
    args.sample_branch_length_script_output_filename = os.path.join(args.output_directory, args.sample_branch_length_script_output_filename)
    args.inference_script_output_filename = os.path.join(args.output_directory, args.inference_script_output_filename)
    case2_script_ancient_samples_file_path = None

    # Actual Computation
    run_step(
        p_initial=args.initial_allele_freq,
        s=args.selection_coefficient,
        n=args.effective_population_size,
        output_file_path=args.step_script_output_file_path,
        ton=args.ton,
        toff=args.toff
    )
    run_mssel(
        nchroms=args.nchroms,
        nreps=args.nreps,
        nder=args.nder,
        nanc=args.nanc,
        path_to_trajfile=args.step_script_output_file_path,
        genome_length=args.genome_length,
        sel_spot=args.sel_spot,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        output_file=args.massel_output,
    )
    convert_txt_to_haps_and_sample(
        r_script_path=args.path_to_converter_script,
        input_txt_file_path=args.massel_output,
        output_file_name=args.converted_filename,
    )
    run_relate(
        relate_binary_path=args.path_to_relate_bin,
        mode=args.relate_mode,
        haps_file_path=str(Path(args.converted_filename).with_suffix(".haps")),
        sample_file_path=str(Path(args.converted_filename).with_suffix(".sample")),
        map_file_path=args.relate_map_file_path,
        effective_population_size=args.effective_population_size*2,
        mutation_rate=args.mutation_rate,
        output_file_path=args.relate_output_filename,
    )
    mut_file_path = os.path.join(args.output_directory, str(Path(args.relate_output_filename).with_suffix('.mut')))
    anc_file_path = os.path.join(args.output_directory, str(Path(args.relate_output_filename).with_suffix('.anc')))

    if os.path.exists(mut_file_path):
        os.remove(mut_file_path)

    if os.path.exists(anc_file_path):
        os.remove(anc_file_path)

    shutil.move(str(Path(args.relate_output_filename).with_suffix('.mut')), args.output_directory)
    shutil.move(str(Path(args.relate_output_filename).with_suffix('.anc')), args.output_directory)

    args.relate_output_filename = os.path.join(args.output_directory, args.relate_output_filename)

    if args.create_ancient_samples:
        case2_script_ancient_samples_file_path = os.path.join(args.output_directory, "ancientSamples.txt")
        run_case2(
            p_initial=args.initial_allele_freq,
            s=args.selection_coefficient,
            n=args.effective_population_size,
            ton=args.ton,
            toff=args.toff,
            ancient_sample_generation_gap=args.case2_script_ancient_samples_generation_gap,
            number_of_ancient_samples=args.case2_script_number_of_ancient_samples,
            output_file_path=case2_script_ancient_samples_file_path
        )

    if args.inference_script_coalescence_times_file_path is not None:
        args.inference_script_coalescence_times_file_path = os.path.join(args.output_directory, args.inference_script_coalescence_times_file_path)

        run_sample_branch_length(
            script_path=args.path_to_sample_branch_length_script,
            input_file_name=args.relate_output_filename,
            mutation_rate=args.mutation_rate,
            coal_file_path=args.sample_branch_length_coal_file_path,
            _format=args.sample_branch_length_script_format,
            output_file_name=args.sample_branch_length_script_output_filename,
            n_samples=args.sample_branch_length_script_n_samples,
        )

    run_inference(
        coalescence_times=args.inference_script_coalescence_times_file_path,
        ancient_samples_file_path=case2_script_ancient_samples_file_path,
        inference_output_filename=args.inference_script_output_filename,
    )
    plot(
        input_file_path=args.inference_script_output_filename,
        output_file_path=os.path.join(args.output_directory, "plot")
    )

if __name__ == "__main__":
    main()
