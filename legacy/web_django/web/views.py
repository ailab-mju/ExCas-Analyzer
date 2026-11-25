from django.shortcuts import render
import os
from django.conf import settings
from django.http import JsonResponse
import json
from django.views.decorators.csrf import csrf_exempt
def index(request):
	return render(request, 'web/index.html')

@csrf_exempt
def sam_upload(request):
	upload_dir=os.path.join(settings.BASE_DIR,"web","static","result_dir")+os.sep
	os.makedirs(upload_dir, exist_ok=True)
	if request.is_ajax:
		sam_list=[]
		sam_name_list=['ctr_only','ctr_com','mut_only','mut_com']
		for sam_name in sam_name_list:
			if request.POST.get(sam_name,None) is not None:
				sam_list.append(json.loads(request.POST.get(sam_name,None)))
			else:
				sam_list.append([])
		for each_sam,sam_name in zip(sam_list,sam_name_list):
			for i in range(0,len(each_sam)):
				each_sam[i] = each_sam[i].split("\t")
			each_sam.sort(key=lambda x : int(x[3]))
			with open(upload_dir+sam_name+".sam", "w") as f:
				for row in each_sam:
					f.write("%s\n" % '\t'.join(str(col) for col in row))
			#print(each_sam)
		context = {
			'ctr_only': sam_list[0],
			'ctr_com': sam_list[1],
			'mut_only': sam_list[2],
			'mut_com': sam_list[3],
		}
		return JsonResponse(context)
	return JsonResponse({'msg':'error'},status=400)
# Create your views here.
