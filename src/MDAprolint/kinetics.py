import pandas as pd

def calculate_residence_times(interactions_per_frame, frame_time):
    """
    Calculates the interaction start time, end time, and continuous duration.
    
    Parameters:
    - interactions_per_frame: List of dicts (output from interactions.py).
    - frame_time: Time elapsed between consecutive analyzed frames in picoseconds (ps).
    
    Returns:
    - Pandas DataFrame of all continuous binding events.
    """
    active_interactions = {} 
    all_events = []          
    
    print("Calculating residence and cumulative times...")
    
    for step_index, frame_data in enumerate(interactions_per_frame):
        current_lipids = frame_data['lipids']
        current_time = step_index * frame_time
        
        for lipid in list(active_interactions.keys()):
            if lipid not in current_lipids:
                start_time = active_interactions.pop(lipid)
                end_time = current_time
                duration = end_time - start_time
                
                all_events.append({
                    'lipid_id': lipid,
                    'start_time': start_time,
                    'end_time': end_time,
                    'duration': duration
                })
                
        for lipid in current_lipids:
            if lipid not in active_interactions:
                active_interactions[lipid] = current_time
                
    final_step_index = len(interactions_per_frame) - 1
    final_time = final_step_index * frame_time
    
    for lipid, start_time in active_interactions.items():
        end_time = final_time
        duration = end_time - start_time
        
        all_events.append({
            'lipid_id': lipid,
            'start_time': start_time,
            'end_time': end_time,
            'duration': duration
        })
        
    df = pd.DataFrame(all_events)
    if not df.empty:
        df = df.sort_values(by=['lipid_id', 'start_time']).reset_index(drop=True)
        
    return df
