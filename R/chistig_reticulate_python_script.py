import pandas as pd 
import random 
import numpy as np
import re


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def extract_numbers_and_convert(text):
    numbers_only = re.sub(r'\D', '', text)  # Remove non-numeric characters
    return int(numbers_only)  # Convert to an integer


def obtain_egos(efile):
	egosdf = pd.read_csv(efile, index_col=0)
	return egosdf

def obtain_egos_list(efile):
	egosdf = pd.read_csv(efile, index_col=0)
	return list(egosdf.egoid)


def egos2apps(efile):
	egosdf = pd.read_csv(efile)
	egosdf = egosdf[['egoid','appslist_rel', 'appslist_norel']]
	#
	egosdf['appslist_rel'].fillna(egosdf['egoid'], inplace=True)
	egosdf['appslist_norel'].fillna(egosdf['egoid'], inplace=True)
	return egosdf

def create_test_rel_status(efile):
	elist = obtain_egos_list(efile)
	#
	rel_status_list = []
	for i in range(len(elist)):
		rel_status_list.append(random.randint(0,1))
	df = pd.DataFrame({'egoid':elist, 'rel_status':rel_status_list})
	return df

 
def categorize_apps(egofile, apptypefile):
	df = pd.read_csv(apptypefile)
	df.rename(columns={'apptype':'appcategory'}, inplace=True)
	df['apptype'] = np.where(df['appcategory'].isin(['hookup-datingapp','classifiedandescort']), 'dating', 'nondating')
	#
	apps_dating = df.loc[df['appcategory'].isin(['hookup-datingapp','classifiedandescort']), 'appid'].tolist()
	apps_nondating = df.loc[~df['appcategory'].isin(['hookup-datingapp','classifiedandescort']), 'appid'].tolist()
	#
	segos2apps = pd.read_csv(egofile)
	segos2apps = segos2apps[['egoid','appslist_rel', 'appslist_norel']]
	segos2apps['appslist_rel'].fillna(segos2apps['egoid'], inplace=True)
	segos2apps['appslist_norel'].fillna(segos2apps['egoid'], inplace=True)
	#
	segos2apps['appslist_rel_dating'] = segos2apps['egoid']
	segos2apps['appslist_rel_nondating'] = segos2apps['egoid']
	segos2apps['appslist_norel_dating'] = segos2apps['egoid']
	segos2apps['appslist_norel_nondating'] = segos2apps['egoid']
	#
	for index, row in segos2apps.iterrows():
		appsrel = row['appslist_rel'].split("|")
		appsnorel = row['appslist_norel'].split("|")
		_egoappsdating_rel = []
		_egoappsnondating_rel = []
		_egoappsdating_norel = []
		_egoappsnondating_norel = []
		#
		if (len(appsrel) == 1) and (appsrel[0] == row['egoid']):
			pass
		else:
			_egoappsdating_rel = [thisapp for thisapp in appsrel if thisapp in apps_dating]
			_egoappsnondating_rel = [thisapp for thisapp in appsrel if thisapp in apps_nondating]
			if _egoappsdating_rel:
				segos2apps.loc[index,'appslist_rel_dating'] = '|'.join(_egoappsdating_rel)
			if _egoappsnondating_rel:
				segos2apps.loc[index,'appslist_rel_nondating'] = '|'.join(_egoappsnondating_rel)
		#
		#
		#
		if (len(appsnorel) == 1) and (appsnorel[0] == row['egoid']):
			pass
		else:
			_egoappsdating_norel = [thisapp for thisapp in appsnorel if thisapp in apps_dating]
			_egoappsnondating_norel = [thisapp for thisapp in appsnorel if thisapp in apps_nondating]
			if _egoappsdating_norel:
				segos2apps.loc[index,'appslist_norel_dating'] = '|'.join(_egoappsdating_norel)
			if _egoappsnondating_norel:
				segos2apps.loc[index,'appslist_norel_nondating'] = '|'.join(_egoappsnondating_norel)
	return segos2apps


def use_apps(relstatus_df, egos2apps_df):
	# check if the egoid is formatted as sXXXXX or is just the integer
	if relstatus_df['egoid'].str.isnumeric().all():
		relstatus_df['egoid'] = 's' + relstatus_df['egoid'].astype(str).str.zfill(5)
	df = pd.merge(egos2apps_df, relstatus_df, on='egoid')
	df['rel_status'] = df['rel_status'].astype(int)
	#
	df['apps_all'] = np.where(df['rel_status'] == 0, df['appslist_norel'], df['appslist_rel'])
	df['apps_dating'] = np.where(df['rel_status'] == 0, df['appslist_norel_dating'], df['appslist_rel_dating'])
	df['apps_nondating'] = np.where(df['rel_status'] == 0, df['appslist_norel_nondating'], df['appslist_rel_nondating'])
	#
	df = df[['egoid','apps_all', 'apps_dating', 'apps_nondating']]
	#
	if df['egoid'].str.isnumeric().all():
		pass
	else:
		df['egoid'] = df['egoid'].apply(extract_numbers_and_convert)
	return df


def obtain_egos_adict(afile):
	empop_venues_matrix = pd.read_csv(afile, index_col=0)
	_venues_list = empop_venues_matrix.columns
	_venues_bool = empop_venues_matrix.apply(lambda x: x > 0)
	# convert the venues dataframe into a dictionary map for each 
	empiricalego2venues = {}
	empiricalego2venues['venue_list'] = _venues_bool.apply(lambda x: list(_venues_list[x.values]), axis=1).to_dict()
	#   
	empiricalego2venues['venue_freq'] = {}
	for ego, venuelist in empiricalego2venues['venue_list'].items():
		empiricalego2venues['venue_freq'][ego] = {}
		for venue in venuelist:
			empiricalego2venues['venue_freq'][ego][venue] = empop_venues_matrix.loc[ego,venue]
	return empiricalego2venues



def categorize_venues(venuetype_file):
	venues_df = pd.read_csv(venuetype_file)
	venues_df.rename(columns={"venuetype":"venuecategory"}, inplace=True)
	venues_df['venuetype'] = np.where(venues_df['venuecategory'].isin(['bar-club','bathhouse']), 'dating', 'nondating')
	venues_df.drop(columns=['venuecategory'], inplace=True)
	return venues_df


def attend_venues(attendance_matrix_file, venues_df):
	venues_dating = venues_df.loc[venues_df['venuetype'] == 'dating', 'venueid'].tolist()
	venues_nondating = venues_df.loc[venues_df['venuetype'] == 'nondating', 'venueid'].tolist()
	
	adict = obtain_egos_adict(attendance_matrix_file)
	venue_attendance = {}
	for thisego, thisegodict in adict['venue_freq'].items():
		_venues_attended = []
		for thisvenue, thisattendancefreq in thisegodict.items():
			for i in range(7):
				if random.random() < thisattendancefreq:
					_venues_attended.append(thisvenue)
		_venues_attended = list(set(_venues_attended))
		_venues_attended = natural_sort(_venues_attended)
		#
		_venues_dating = [thisvenue for thisvenue in _venues_attended if thisvenue in venues_dating]
		_venues_nondating = [thisvenue for thisvenue in _venues_attended if thisvenue in venues_nondating]
		#
		venue_attendance[thisego] = {}
		if _venues_attended:
			venue_attendance[thisego]['all'] = '|'.join(_venues_attended)
		else:
			venue_attendance[thisego]['all'] = thisego
		if _venues_dating:
			venue_attendance[thisego]['dating'] = '|'.join(natural_sort(_venues_dating))
		else:
			venue_attendance[thisego]['dating'] = str(thisego)
		if _venues_nondating:
			venue_attendance[thisego]['nondating'] = '|'.join(natural_sort(_venues_nondating))
		else:
			venue_attendance[thisego]['nondating'] = str(thisego)

	egoslist = []
	venuelist = []
	venuelist_dating = []
	venuelist_nondating = []
	for thisego, thisdict in venue_attendance.items():
		egoslist.append(thisego)
		venuelist.append(thisdict['all'])
		venuelist_dating.append(thisdict['dating'])
		venuelist_nondating.append(thisdict['nondating'])

	df = pd.DataFrame({'egoid':egoslist, 'venues_all':venuelist, 'venues_dating':venuelist_dating, 'venues_nondating':venuelist_nondating})
	if df['egoid'].str.isnumeric().all():
		pass
	else:
		df['egoid'] = df['egoid'].apply(extract_numbers_and_convert)
	return df

