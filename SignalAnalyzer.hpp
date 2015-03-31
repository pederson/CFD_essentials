#ifndef _SIGNALANALYZER_H
#define _SIGNALANALYZER_H

#include <iostream>
#include <vector>
#include <string>


class SignalAnalyzer{
public:

	const double & result() const {return m_result.front();};

	void bind_signal(const double * signal, unsigned int length);
	void run_filter();
	void run_FFT();

private:
	const double * m_signal;
	unsigned int m_signal_length;

	std::vector<double> m_result;
};

#endif