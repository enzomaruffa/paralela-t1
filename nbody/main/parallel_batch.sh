export TIMEFORMAT=%R
INPUT_SIZES=(100 1000 5000 10000 15000 25000 40000)
TIMESTEPS=(10 50)
THREADS=(1 2 4 8 12)
ITERATIONS=30

mkdir -p .temp
mkdir -p data

printf "particle_count,timesteps,threads,sample,time\n" > data/times_results.txt

for INPUT in ${INPUT_SIZES[*]}; do
    for TIMESTEP in ${TIMESTEPS[*]}; do
        for THREAD_COUNT in ${THREADS[*]}; do
            export OMP_NUM_THREADS=$THREAD_COUNT
            for ITERATION in `seq $ITERATIONS`; do
                echo "Running ($INPUT, $TIMESTEP, $THREAD_COUNT) iteration $ITERATION"
                
                printf "${INPUT}\n${TIMESTEP}\n" > .temp/parallel_batch_input.in

                RUN_TIME=$({ time ./nbody < .temp/parallel_batch_input.in >.temp/temp_output.txt; } 2>&1)

                diff .temp/temp_output.txt ../original/outputs/${INPUT}_${TIMESTEP}.txt > .temp/temp_diff.txt

                if [ -s .temp/temp_diff.txt ]; then
                    # The file is not-empty. Add to errors.
                    mkdir -p errors
                    echo "  - Error!"
                    printf "Error for input ${INPUT}/${TIMESTEP} with ${THREAD_COUNT} threads in iteration ${ITERATION}\n" >> errors/errors.log
                    # mv .temp/temp_output.txt errors/${INPUT}_${TIMESTEP}_$THREAD_COUNT.out
                    # mv .temp/temp_diff.txt errors/${INPUT}_${TIMESTEP}_$THREAD_COUNT.err
                fi

                printf "${INPUT},${TIMESTEP},${THREAD_COUNT},${ITERATION},${RUN_TIME}\n" >> data/times_results.txt
            done
        done
    done
done

rm .temp/parallel_batch_input.in
rm .temp/temp_output.txt
rm .temp/temp_diff.txt
rm -rf .temp