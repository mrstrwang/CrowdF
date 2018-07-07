from datetime import datetime, tzinfo, timedelta

from django.shortcuts import render, redirect
from django.views.generic.base import View
from django.http import HttpResponse
from django.urls import reverse

import pytz
from operations.models import UserInterest, EmailVerifyCode
from projects.models import Project
from utils.SendVerifyCode import send_email_code
from users.models import UserProfile


class InterestManageView(View):
    def post(self, request, project_id):
        has_interest = UserInterest.objects.filter(user_id=request.user.id, project_id=project_id)
        if not has_interest:
            UserInterest.objects.create(user_id=request.user.id, project_id=project_id)
            self.update(project_id, 1)
            return HttpResponse('{"res":0}', content_type='application/json')
        else:
            has_interest.delete()
            self.update(project_id, -1)
            return HttpResponse('{"res":1}', content_type='application/json')

    def update(self, project_id, num):
        project = Project.objects.get(id=project_id)
        project.interest_nums += num
        project.save()


class AccountVerifyView(View):
    def get(self, request, step_id):
        step_id_relationship = {0: 'member.html',
                                1: 'apply.html',
                                2: 'apply-1.html',
                                3: 'apply-2.html',
                                4: 'apply-3.html'}
        if step_id == 0 and request.user.account_status == 'non_choice':
            step_id_relationship[step_id] = 'accttype.html'
        if step_id == 1:
            request.session['verify_type'] = request.GET.get('verify_type', '')
        return render(request, step_id_relationship[step_id])

    def post(self, request, step_id):
        if step_id == 2:
            for item in ['realname', 'idcard_num', 'mobile']:
                request.session[item] = request.POST.get(item, '')
            # realname = request.POST.get('realname', '')
            # idcard_num = request.POST.get('idcard_num', '')
            # mobile = request.POST.get('mobile', '')
            # request.session['realname'] = realname
            # request.session['idcard_num'] = idcard_num
            # request.session['mobile'] = mobile
            return render(request, 'apply-1.html')

        if step_id == 3:
            # 保存照片
            return render(request, 'apply-2.html')

        if step_id == 4:
            # 提交邮箱，进行确认
            if request.POST.get('sendmail'):
                email = request.POST.get('email', '')
                request.session['email'] = email
                send_email_code(email, send_type='verify_id')
                return render(request, 'apply-3.html', {'email':email})
            else:
                verify_code = request.POST.get('verify_code', '')
                email = request.session['email']
                if verify_code:
                    vaild = EmailVerifyCode.objects.filter(verify_code=verify_code, send_type='verify_id', email=email)
                    if vaild:
                        user = UserProfile.objects.get(id=request.user.id)
                        user.realname = request.session['realname']
                        user.id_card_num = request.session['idcard_num']
                        user.mobile = request.session['mobile']
                        user.email = email
                        user.account_status = request.session['verify_type']
                        user.account_type = 'company'
                        user.save()
                        return render(request, 'member.html')
                else:
                    return render(request, 'apply-3.html', {'error':'验证码错误','email':email})


class RefreshProjectDataView(View):
    def get(self, request):
        all_project = Project.objects.all()
        for project in all_project:
            target_day = project.create_time + timedelta(days=project.deadline)
            today = pytz.utc.localize(datetime.now())
            project.left_days = (target_day - today).days + 1
            project.due_date = target_day
            if project.support_fund > project.target_fund:
                project.status = 'success'
            elif project.left_days < 1:
                project.left_days = 0
                project.status = 'failed'
            else:
                project.status = 'funding'
            project.save()
        return redirect(reverse('index'))
