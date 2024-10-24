import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import Angle
from astropy import units as u
from astropy.coordinates import SkyCoord,FK5,get_sun,EarthLocation
from astropy.time import Time,TimeDelta
import argparse
import os
import psrqpy
import lofarantpos.db
import astroplan
import re

defaultScope = 'IE613'
pulsarNameRE = re.compile(r"[BJ]\d{4}[-+]\d{2,4}", re.IGNORECASE)
db = psrqpy.QueryATNF()

class Observable(object):

	def __init__(self, sourceName, obsStartTime, obsEndTime = None, obsLength = None, ra = None, dec = None, stationName = defaultScope, warningAlt = 40., cutoffAlt = 0.,durn=None):


		self.stationLocationXYZ = lofarantpos.db.LofarAntennaDatabase().phase_centres[stationName + 'HBA']
		self.stationLocation = EarthLocation.from_geocentric(self.stationLocationXYZ[0], self.stationLocationXYZ[1], self.stationLocationXYZ[2], unit = u.m)
		self.obsWarningAltitude = warningAlt * u.deg
		self.obsCutoffAltitude = cutoffAlt * u.deg
		self.durn = durn

		self.obsAPObj = None
		self.obsTarget = None

		if isinstance(ra, type(None)):
			self.setNamedTarget(sourceName)
		else:
			self.setTarget(sourceName, ra, dec)

		self.setObsTime(obsStartTime, obsEndTime, obsLength * u.hour)





	def setNamedTarget(self, sourceName):
		pulsar_name = sourceName
		if not isinstance(pulsarNameRE.match(sourceName), type(None)):
			sourceName = f"{sourceName}_PSR"
		#self.obsAPObj = SkyCoord.from_name(sourceName)
		self.obsAPObj = SkyCoord(db[pulsar_name]['RAJ'].value[0],db[pulsar_name]['DECJ'].value[0], unit = 'hourangle,degree')
		print(f"Set target to {sourceName} at ({self.obsAPObj.to_string('hmsdms')})")
		self.obsTarget = sourceName


	def setTarget(self, sourceName, ra, dec, units = (u.hourangle, u.deg)):
		self.obsAPObj = SkyCoord(ra = ra, dec = dec, unit = units, frame = 'icrs')
		#print(f"Set target to {sourceName} at ({self.obsAPObj.to_string('hmsdms')})")
		self.obsTarget = sourceName


	def setObsTime(self, obsStartTime = None, obsEndTime = None,obsLength = 1 * u.day):
		if not isinstance(obsStartTime, type(None)):
			self.obsStartTime = Time(obsStartTime, format = 'isot')
		else:
			self.obsStartTime = Time.now()
				

		if not isinstance(obsEndTime, type(None)):
			self.obsEndTime = Time(obsEndTime, format = 'isot')
		else:
			self.obsEndTime = self.obsStartTime + obsLength


		#print(f"Set start / end of window to {self.obsStartTime} - {self.obsEndTime}")


	def checkObservationTarget(self, savePlot = False, outputPath = "./"):
		"""Check if the target falls out of view or below configured cutoffs during the observations
		
		Returns:
		    bool: Success
		
		Raises:
		    RuntimeError: Observation not properly configured
		"""
		# name resolved, alt/az
		if self.obsAPObj is None or self.obsStartTime is None or self.obsEndTime is None:
			raise RuntimeError(f"Observation Improperly Defined. Start/End Time: {self.obsStartTime}/{self.obsEndTime}, APObj: {self.obsAPObj}.")

		aplanObs = astroplan.Observer(location = self.stationLocation)
		warningAlt = astroplan.AltitudeConstraint(min = self.obsWarningAltitude)
		cutoffAlt = astroplan.AltitudeConstraint(min = self.obsCutoffAltitude)

		altData = []
		currTime = self.obsStartTime
		timeDelta = TimeDelta(15 * 60, format='sec')
		duration = TimeDelta(float(self.durn)*60, format='sec')
		while currTime < self.obsEndTime:
			alt = aplanObs.altaz(currTime, self.obsAPObj).alt.deg
			altData.append((currTime.plot_date, alt, float(alt)))
			currTime += timeDelta

		altData = np.vstack(altData)
		plt.figure(figsize = (24,12))
		plt.plot_date(altData[:, 0], altData[:, 1], label = f'Altitude of {self.obsTarget} (deg)')
		yLim = plt.ylim()
		xLim = plt.xlim()

		warningTransition = np.where(np.diff(np.sign(altData[:, 2] - self.obsWarningAltitude.value)))
		cutoffTransition = np.where(np.diff(np.sign(altData[:, 2] - self.obsCutoffAltitude.value)))
		fracTime = (altData[-1, 0] - altData[0, 0]) / 100
		if (warningTransition[0].size > 0):

			for idx in np.nditer(warningTransition):
				plt.vlines(altData[idx, 0], yLim[0], yLim[1], color = 'y', alpha = 0.4)
				time = altData[idx, 0] + 0.33 * fracTime
				if altData[idx, 2] - self.obsWarningAltitude.value > 0:
					time -= 1.27 * fracTime
				plt.text(time,  np.mean(yLim), "(W) " + str(self.obsStartTime + idx * timeDelta), rotation = 90, verticalalignment = 'center', color = 'gray')


			if (cutoffTransition[0].size > 0):
				for idx in np.nditer(cutoffTransition):
					plt.vlines(altData[idx, 0], yLim[0], yLim[1], color = 'r', alpha = 0.66)
					time = altData[idx, 0] + 0.33 * fracTime
					if altData[idx, 2] - self.obsCutoffAltitude.value< 0:
						time -= 1.27 * fracTime
					plt.text(time, np.mean(yLim), "(C) " + str(self.obsStartTime + idx * timeDelta), rotation = 90, verticalalignment = 'center', color = 'black')

		peakAlt = np.argmax(altData[:, 2])
		plt.vlines(altData[peakAlt, 0], yLim[0], yLim[1], color = 'g', alpha = 0.4)
		time = altData[peakAlt, 0] + 0.33 * fracTime
		plt.text(time,  np.mean(yLim), "(P) " + str(self.obsStartTime + peakAlt * timeDelta), rotation = 90, verticalalignment = 'center', color = 'gray')

		#print(str(self.obsStartTime + peakAlt * timeDelta - duration/2)+" - "+str(self.obsStartTime + peakAlt * timeDelta + duration/2)+" -\t"+str(self.obsTarget))
		st_time=self.obsStartTime + peakAlt * timeDelta - duration/2
		ed_time=self.obsStartTime + peakAlt * timeDelta + duration/2
		st_time=st_time.strftime('%Y-%m-%dT%H:%M')
		ed_time=ed_time.strftime('%Y-%m-%dT%H:%M')
		print(f"{st_time} - {ed_time} : {self.obsTarget} [{self.obsAPObj.ra.rad},{self.obsAPObj.dec.rad}, \'J2000\']")
        

		titleText = f"{self.obsTarget}@{self.obsStartTime}-{self.obsEndTime}\nMin {np.min(altData[:, 1]):.2f}deg, Max {np.max(altData[:, 1]):.2f}deg"
		minAlt = np.min(altData[:, 1])
		if minAlt < self.obsCutoffAltitude.value:
			plt.fill_between(xLim, yLim[0], self.obsCutoffAltitude.value, color = 'r', alpha = 0.33, label = f"Cutoff Altitude ({self.obsCutoffAltitude.value}$^\circ$)")				
			titleText += " CUFOFF"
		if minAlt < self.obsWarningAltitude.value:
			plt.fill_between(xLim, max(minAlt, self.obsCutoffAltitude.value), self.obsWarningAltitude.value, color = 'y', alpha = 0.2, label = f"Warning Altitude ({self.obsWarningAltitude.value}$^\circ$)")
			titleText += " WARNING"

		plt.ylabel("Altitude (deg)")
		plt.xlabel("Time (UTC)")
		plt.title(titleText)
		plt.xlim(xLim)
		plt.ylim(yLim)
		plt.legend()
		if savePlot:
			plt.savefig(os.path.join(outputPath, f"altData_{self.obsTarget}_{self.obsStartTime}_{self.obsEndTime}.png"))
		#plt.show()



		alwaysObservable = np.array(astroplan.is_always_observable(warningAlt, aplanObs, self.obsAPObj, time_range = [self.obsStartTime, self.obsEndTime]))
		if alwaysObservable.all():
			return altData
		else:
			print("")
			visibleFrac = np.sum(altData[:, 1] < self.obsWarningAltitude.value) * 100  / altData.shape[0]
			#print(f"Altitude falls below warning altitude {self.obsWarningAltitude} for {visibleFrac:3.1f}% of the observation window.")

			invisibleFrac = np.sum(altData[:, 1] < self.obsCutoffAltitude.value) * 100  / altData.shape[0]
			if invisibleFrac > 0:
				print("")#print(f"Target falls below cuttoff altitude {self.obsCutoffAltitude} for {invisibleFrac:3.1f}% of the observation window.")

		return altData





if __name__== '__main__':
	parser = argparse.ArgumentParser(description = "Plot the altitude for a source over time.")
	parser.add_argument('--source', '-s', dest = 'src', required = True, type = str,nargs='+', help = "Name of the source(s)")
	parser.add_argument('--ra', '-r', dest = 'ra', default = None,  type = str, help = "Right ascension of a source in J2000 /  HH:MM:SS.ss format.")
	parser.add_argument('--dec', '-d', dest = 'dec', default = None,  type = str, help = "Right ascension of a source in J2000 /  HH:MM:SS.ss format.")
	parser.add_argument('--start', '-t', dest = 'start', default = None, type = str, help = "ISOT start time in the YYYY-MM-SSTHH:MM:SS format.")
	parser.add_argument('--end', '-e', dest = 'end', default = None, type = str, help = "ISOT end time in the YYYY-MM-SSTHH:MM:SS format.")
	parser.add_argument('--length', '-l', dest = 'length', type = float, default = 24, help = "Amount of time to scan from start time in hours from start (replaces --end)")
	parser.add_argument('--tele', dest = 'telescope', default = defaultScope, help = "LOFAR station ID, e.g. IE613")
	parser.add_argument('--warn', dest = 'warn', default = 40, type = float, help = "Altitude to apply warning layer at (deg).")
	parser.add_argument('--cutoff', dest = 'cut', default = 0, type = float, help = "Altitude to apply cutoff layer at (deg).")
	parser.add_argument('--plot_only', '-p', default = False, dest = 'plot', action = 'store_true', help = "Only plot the data, don't write figure to disk.")
	parser.add_argument('--output', '-o', dest = 'output', default = './', help = "Location to save plots (if not using --plot_only)")
	parser.add_argument('--durn', '-durn', dest = 'duration', default = None, help = "Required Observation duration of the source")    

	args = parser.parse_args()
	for src in args.src:
		obs = Observable(sourceName=src, obsStartTime=args.start, obsEndTime=args.end, obsLength=args.length, ra=args.ra, dec=args.dec, stationName=defaultScope, warningAlt=args.warn, cutoffAlt=args.cut, durn=args.duration)
		obs.checkObservationTarget(not args.plot, args.output)








