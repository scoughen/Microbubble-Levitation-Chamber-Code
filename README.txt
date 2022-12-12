The following is a summary of the purpose/use of each of the programs contained in this repository.

LabVIEW Code: 
	Custom Oscope Functions:
		Oscop 2 Channel Phase Diff.vi - This vi reads and outputs the phase difference between the two oscilloscope channels.
		
		Oscop Phase Read Averaging.vi - This vi uses "Oscop 2 Channel Phase Diff.vi" to take a series of phase readings and then average them.
		
		Oscop Read Averaging.vi - This vi uses "Oscop Read.vi" to take a series of amplitude readings and then average them.
		
		Oscop Read.vi - This vi reads and outputs the amplitude (or some other values) from a specified oscilloscope channel.
		
		Phase Diff Setup.vi - This vi performs the necessary setup for using the oscilloscope phase read vi's.
		
		
	Marlin Instrument Driver:
		This is the driver used for controlling the Ender 3 3D printer.
		
		
	Siglent SDS 1000 2000 Series:
		This is the driver used for interfacing with the Siglent 1202X-E oscilloscope.
		
		
	Read Waveform (Single) With Data Output.vi - This vi reads the entire measured waveform on the oscilloscope for the specified channel and outputs the resulting data to a .csv file.
	
	Scanning Motion One Direction for Y Multiple Readings per Point With Phase Info.vi - This vi uses the Ender 3 and the oscilloscope to perform an acoustic scan for the specified parameters. It returns both amplitude and phase data.
	
	Scanning Motion One Direction for Y Multiple Readings per Point.vi - This vi uses the Ender 3 and the oscilloscope to perform an acoustic scan for the specified parameters. It returns only amplitude data.
	
	

FFT.m - This script takes the time series of a signal read from an oscilloscope and calculates/plots the FFT.

HorizontalScanPlottingPlusPhase.m - This script takes the data (both amplitude and phase) from a horizontal (xy-plane) acoustic scan and plots it.

ImageProcessing.m - This script takes an image of a (spherical) bubble and extracts the position and radius.

NonsphericalOscillationAnalysis_MLC.m - This script takes an image of a bubble, extracts the r(theta) radial description of the bubble, and calculates the nonspherical shape modes of the bubble.

ScanAndTheoryComparison.m - This script plots the scan of the vertical acoustic field of a plane transducer, calculates the theoretical acoustic field for a plane transducer of the given specs and driving conditions (these should be set to match thoes of the scanned transducer), and compares them.

TankModesCalculations.m - This scripts calculates the modes of a rectangular tank for the given dimensions. The resulting modes and plotted.

TheoreticalTransducerRadiationCalcs_FinalVersion.m - This script calculates the theoretical acoustic field of a plane transducer for the given specs and driving conditions.

VerticalScanPlottingPlusPhase.m - This script takes the data (both amplitude and phase) from a veritical acoustic scan in the xz-plane and plots it.

VerticalScanPlottingPlusPhaseOrthogonal.m - This script takes the data (both amplitude and phase) from a vertical acoustic scan in the yz-plane and plots it.

VideoProcessing.m - This script takes a video of a bubble and extracts the position and radius for each frame and plots the resulting radius time series.
