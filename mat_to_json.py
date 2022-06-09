import cobra

model = cobra.io.load_matlab_model('model.mat')
cobra.io.save_json_model(model, 'model.json')