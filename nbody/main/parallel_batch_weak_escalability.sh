export TIMEFORMAT=%R
INPUT_SIZES=(1000 2000 4000 8000 12000)
TIMESTEPS=(50)
THREADS=(1 2 4 8 12)
ITERATIONS=30

mkdir -p .temp
mkdir -p data

printf "particle_count,timesteps,threads,sample,time\n" > data/weak_escalability_times_results.txt

for INPUT in ${INPUT_SIZES[*]}; do
    for TIMESTEP in ${TIMESTEPS[*]}; do
        for THREAD_COUNT in ${THREADS[*]}; do
            export OMP_NUM_THREADS=$THREAD_COUNT
            for ITERATION in `seq $ITERATIONS`; do
                echo "Running ($INPUT, $TIMESTEP, $THREAD_COUNT) iteration $ITERATION"
                
                printf "${INPUT}\n${TIMESTEP}\n" > .temp/parallel_batch_input.in

                RUN_TIME=$({ time ./nbody < .temp/parallel_batch_input.in >.temp/temp_output.txt; } 2>&1)

                printf "${INPUT},${TIMESTEP},${THREAD_COUNT},${ITERATION},${RUN_TIME}\n" >> data/weak_escalability_times_results.txt
            done
        done
    done
done

rm .temp/parallel_batch_input.in
rm .temp/temp_output.txt
rm .temp/temp_diff.txt
rm -rf .temp