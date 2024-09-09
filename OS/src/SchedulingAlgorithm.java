// Tala Atallah -- 1202575
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Comparator;

public class SchedulingAlgorithm {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int previousArrival = 0; // this variable is used to assign some kind of order for arrival time
		Random random = new Random();
		int iterations = 100000; // initialize the iterations

		// i created lists for each algorithm ,one for awt and one for tat , to store
		// the values
		List<Double> awtListFCFS = new ArrayList<>();
		List<Double> tatListFCFS = new ArrayList<>();
		List<Double> awtListSRTF = new ArrayList<>();
		List<Double> tatListSRTF = new ArrayList<>();
		List<Double> awtListRR = new ArrayList<>();
		List<Double> tatListRR = new ArrayList<>();
		List<Double> awtListMFQ = new ArrayList<>();
		List<Double> tatListMFQ = new ArrayList<>();

		for (int i = 1; i <= iterations; i++) {
			List<int[]> processes = new ArrayList<>(); // Create an array list of processes that contains arrays of
														// integer
			int arrivalTime;
			int CPUburst;
			for (int j = 0; j < 8; j++) { // for loop for creating 8 process
				int[] process = new int[3];
				process[0] = j; // index 0 is assigned as an id for the process
				arrivalTime = (int) (previousArrival * Math.random()) + 2 * j; // the function i created to assign
																				// arrival time in some random order
				process[1] = arrivalTime; // index 1 is the arrival time
				CPUburst = random.nextInt(96) + 5; // randomly generate the cpu time between 5 and 100
				process[2] = CPUburst; // assign index 2 to the cpu burst time
				processes.add(process); // then add the array (process) to the arraylist that contains all the 8
										// processes
				previousArrival = arrivalTime;
			}

			processes.sort(Comparator.comparingInt(arr -> arr[1])); // sort the processes based on their arrival time

			if (i == 100 || i == 1000 || i == 10000 || i == 100000) { // check if the iterations are equal to
				// 100,1000,10000,100000 to print the results

				double[] resultsFCFS = FCFS(processes); // call the function FCFS and store there returned values
				// add the tat and awt to their lists
				awtListFCFS.add(resultsFCFS[0]);
				tatListFCFS.add(resultsFCFS[1]);

				double[] resultsSRTF = SRTF(processes); // call the function SRTF and store there returned values
				// add the tat and awt to their lists
				awtListSRTF.add(resultsSRTF[0]);
				tatListSRTF.add(resultsSRTF[1]);

				double[] resultsRR = RR(processes); // call the function RR and store there returned values
				// add the tat and awt to their lists
				awtListRR.add(resultsRR[0]);
				tatListRR.add(resultsRR[1]);

				double[] resultsMFQ = MFQ(processes); // call the function MFQ and store there returned values
				// add the tat and awt to their lists
				awtListMFQ.add(resultsMFQ[0]);
				tatListMFQ.add(resultsMFQ[1]);

			}

		}
		// call the function print table for the 4 algorithms
		printTable("FCFS", awtListFCFS, tatListFCFS);
		printTable("SRTF", awtListSRTF, tatListSRTF);
		printTable(" RR ", awtListRR, tatListRR);
		printTable("MLFQ", awtListMFQ, tatListMFQ);

	}

	private static double[] FCFS(List<int[]> processes) {

		int currTime = 0;
		int finishTime;
		int turnaroundTime = 0;
		int waitingTime = 0;

		for (int i = 0; i < processes.size(); i++) { // loop through each process
			int[] process = processes.get(i);
			int arrivalTime = process[1]; // assign arrival time to index 1
			int CPUBurst = process[2]; // assign cpu burst time to index 2

			finishTime = Math.max(currTime, arrivalTime) + CPUBurst;
			turnaroundTime += finishTime - arrivalTime; // calculate the accumulative TAT for all processes by
														// subtracting finish time from arrival time
			waitingTime += finishTime - arrivalTime - CPUBurst; // calculate the accumulative WT for all processes by
																// subtracting finish time from arrival time and cputime
			currTime = finishTime;

		}
		double AWT = 0;
		double TAT = 0;
		// after the loop finishes calculate the average turnaround time and waiting
		// time by dividing the accumulative values by the number of processes
		AWT = (double) waitingTime / 8;
		TAT = (double) turnaroundTime / 8;
		return new double[] { AWT, TAT }; // return the values of average waiting time and average turn around time to
											// the main.

	}

	private static double[] SRTF(List<int[]> processes) {
		int[] tempBurst = new int[8];
		int[] arrivalTime = new int[8];
		int[] burstTime = new int[8];
		int[] completedTime = new int[8];
		int[] waitingTime = new int[8];

		int sumwaitingTime = 0;
		int currTime = 0;
		int turnaroundTime = 0;
		int completedProcess = 0;
		for (int i = 0; i < processes.size(); i++) { // loop through each process
			int[] process = processes.get(i);
			arrivalTime[i] = process[1]; // assign arrival time to index 1
			burstTime[i] = process[2]; // assign cpu burst time to index 2
			tempBurst[i] = burstTime[i]; // assign a copy of cpu burst time
		}
		while (completedProcess < 8) { // the process will
			int minProcessIndex = -1; // initially assign the index number that has the least cpu time to -1
			int minCPUTime = Integer.MAX_VALUE; // and a value for the actual minimum cpu time
			for (int i = 0; i < 8; i++) { // loop through the 8 processes

				if (tempBurst[i] > 0 && arrivalTime[i] <= currTime && tempBurst[i] < minCPUTime) {
					minCPUTime = tempBurst[i];
					minProcessIndex = i;
				}

			}
			if (minProcessIndex != -1) { // if the minimum index is not -1 
				tempBurst[minProcessIndex]--; // then decrement the burst time 
				if (tempBurst[minProcessIndex] == 0) { // check if its equal to zero this means it has finished

					completedTime[minProcessIndex] = currTime; // assign its completion time to the current time
					// calculate its waiting time , using math.max(0,..) ensures that the waiting time won't be negative
					waitingTime[minProcessIndex] = Math.max(0,completedTime[minProcessIndex] - arrivalTime[minProcessIndex] - burstTime[minProcessIndex]);

					completedProcess++; // increment how many processes are completed 
				}

			}
			currTime++; // then update the current time

		}
		// for loop to do the summation of all waiting times and turn around times
		for (int i = 0; i < processes.size(); i++) {
			int[] process = processes.get(i);
			turnaroundTime += waitingTime[i] + process[2];

			sumwaitingTime += waitingTime[i];
		}
		//calculate the average
		double AWT;
		double TAT;
		AWT = (double) sumwaitingTime / 8;
		TAT = (double) turnaroundTime / 8;
		return new double[] { AWT, TAT };

	}

	private static double[] RR(List<int[]> processes) {
		int quantum = 20; // set the quantum to 20
		int leftProcesses = 8;
		int[] tempBurst = new int[8]; // tempburst is used to store the burst that is been decremented until it
										// reaches zero
		int[] tempArrival = new int[8];
		int[] burstTime = new int[8];
		int currTime = 0;
		int waitingTime = 0;
		int turnaroundTime = 0;
		for (int i = 0; i < processes.size(); i++) {
			int[] process = processes.get(i);
			tempArrival[i] = process[1];
			burstTime[i] = process[2];
			tempBurst[i] = burstTime[i];
		}
		int j = 0; // starting with process zero
		int flag = 0;
		while (leftProcesses != 0) { // while loop : checks if there are remaining processes that are not fully
										// excuted
			if (tempBurst[j] <= quantum && tempBurst[j] > 0) { // checks if the burst of the current process is less
																// than quantum
				currTime += tempBurst[j]; // update the current time
				tempBurst[j] = 0; // set the temp burst to zero , so it won't be negative
				flag = 1; // this flag is to indicate that the process is done
			} else if (tempBurst[j] > 0) { // if the burst is larger than zero
				tempBurst[j] = tempBurst[j] - quantum; // subtract the quantum from the burst time
				currTime += quantum; // update the current time
			}

			if (flag == 1 && tempBurst[j] == 0) { // this if statement is to ensure that the process is done executing
				leftProcesses = leftProcesses - 1; // if the condition is true then we will decrement the number of
													// processes left to be executing

				turnaroundTime = turnaroundTime + currTime - tempArrival[j]; // calculate the turnaround time
				waitingTime = waitingTime + currTime - tempArrival[j] - burstTime[j]; // calculate the waiting time
				flag = 0; // and then set the flag to zero again for the other processes
			}
			if (j == 7) { // this condition checks if we reached to process the last process if so then we
							// should go back to execute the first processes
				j = 0;
			} // update the index to zero ( indicates process one)

			else if (tempArrival[j + 1] <= currTime) { // if the arrival time of the next process is less than the
														// current time
				j++;
			} // then increment the index and go to the next process to
				// be executed

			else {
				j = 0;
			}

		}
		// calculate the average waiting and turn around time by dividing by 8
		double AWT;
		double TAT;
		AWT = (double) waitingTime / 8;
		TAT = (double) turnaroundTime / 8;
		return new double[] { AWT, TAT };

	}

	private static double[] MFQ(List<int[]> processes) {

		int[] tempBurst = new int[8];
		int[] tempArrival = new int[8];
		int[] burstTime = new int[8];
		int currTime = 0;
		int[] completedTime = new int[8];
		int waitingTime = 0;
		int turnaroundTime = 0;

		int[] quantum = { 10, 50 }; // set values of quantums in RR

		// create linked lists ( acts as queues , each with different algorithm
		// implemented in )
		LinkedList<int[]> q1 = new LinkedList<>();
		LinkedList<int[]> q2 = new LinkedList<>();
		LinkedList<int[]> q3 = new LinkedList<>();
		// in this for loop we'll iterate through the array list to access the arrays
		// (processes) inserted in
		for (int i = 0; i < processes.size(); i++) {
			int[] process = processes.get(i);
			tempArrival[i] = process[1];
			burstTime[i] = process[2];
			tempBurst[i] = burstTime[i];
			q1.add(process); // add each process in queue 1

		}

		while (!q1.isEmpty() || !q2.isEmpty() || !q3.isEmpty()) { // the loop will remain executing until all the queues
																	// are empty , this means all the processes were
																	// done
			// the first queue will implement RR with quantum 10
			for (int i = 0; i < q1.size(); i++) {
				int[] currentProcess = q1.removeFirst(); // start by accessing the first process
//				System.out.println("Process " + currentProcess[0] + " running from queue 1 at time " + currTime);
				if (currentProcess[2] <= quantum[0]) { // check if the cpu burst time for the current process is less
														// than the first quantum
					currTime += currentProcess[2]; // if so thwn update the current time
					currentProcess[2] = 0; // and assign the cpu burst time to zero
				} else {
					currTime += quantum[0]; // else update the current time with the quantum
					currentProcess[2] -= quantum[0]; // subtract the cpuburst time with the quantum
					q2.add(currentProcess); // then add the process to the next queue ( it will be implemented in the
											// next queue if it didn't finish execution in the first queue)
				}

				if (currentProcess[2] <= 0) { // check if the cpu burst is equal to zero to calculate its turnaround and
												// waiting time
					completedTime[currentProcess[0]] = currTime; // currentProcess[0] is the id of the process , its
																	// completion time will be updated to the current
																	// time

					// calculate the waiting time by subtracting burst time and arrival time from
					// the completion time
					waitingTime = completedTime[currentProcess[0]] - currentProcess[1] - burstTime[currentProcess[0]];

					turnaroundTime = completedTime[currentProcess[0]] - currentProcess[1]; // calculate turn around time
																							// by subtracting arrival
																							// time from completion time
				}

			}
			// check that all the processes are implemented in the first queue
			if (q1.isEmpty()) {
				// the following will be the same as the first queue
				// RR is implemented but the quantum is different
				for (int i = 0; i < q2.size(); i++) {
					int[] currentProcess = q2.removeFirst();
					if (currentProcess[2] <= quantum[1]) {
						currTime += currentProcess[2];
						currentProcess[2] = 0;
					} else {
						currTime += quantum[1];
						currentProcess[2] -= quantum[1];
						q3.add(currentProcess); // add the processes that are not done executing to the last queue
					}

					if (currentProcess[2] <= 0) {
						completedTime[currentProcess[0]] = currTime;
						waitingTime += completedTime[currentProcess[0]] - currentProcess[1]
								- burstTime[currentProcess[0]];
						turnaroundTime += completedTime[currentProcess[0]] - currentProcess[1];
					}

				}
			}
			// check if the processes in q2 are fully excuted
			if (q2.isEmpty()) {
				// q3 will run fcfs algorithm
				for (int i = 0; i < q3.size(); i++) {
					int[] currentProcess = q3.removeFirst();
					currTime += currentProcess[2]; // it will execute the left cpu burst left
					completedTime[currentProcess[0]] = currTime; // update completion time
					waitingTime += completedTime[currentProcess[0]] - currentProcess[1] - burstTime[currentProcess[0]]; // calculate
																														// waiting
																														// time
					turnaroundTime += completedTime[currentProcess[0]] - currentProcess[1]; // calculate turn around
																							// time
				}
			}
		}
		// average turn around time and waiting time are calculated by dividing by 8 ( # of processes)
		double AWT;
		double TAT = 0;
		AWT = (double) waitingTime / 8;
		TAT = (double) turnaroundTime / 8;

		return new double[] { AWT, TAT };

	}

	// the following is to print as a table , by taking the name of the algorithm and the lists 
	// that stores the awt and tat for each iteration
	public static void printTable(String name, List<Double> awtList, List<Double> tatList) {
		System.out.println("================================= " + name + " ================================= ");
		System.out.println("itr\t\t100\t\t1000\t\t10000\t\t100000");

		if (!awtList.isEmpty()) {
			System.out.print("AWT\t\t" + awtList.get(0));
		}
		if (awtList.size() > 1) {
			System.out.print("\t\t" + awtList.get(1));
		}
		if (awtList.size() > 2) {
			System.out.print("\t\t" + awtList.get(2));
		}
		if (awtList.size() > 3) {
			System.out.print("\t\t" + awtList.get(3));
		}
		System.out.println();

		if (!tatList.isEmpty()) {
			System.out.print("TAT\t\t" + tatList.get(0));
		}
		if (tatList.size() > 1) {
			System.out.print("\t\t" + tatList.get(1));
		}
		if (tatList.size() > 2) {
			System.out.print("\t\t" + tatList.get(2));
		}
		if (tatList.size() > 3) {
			System.out.print("\t\t" + tatList.get(3));
		}
		System.out.println();

	}
}
