from django.shortcuts import render, redirect
from django.http import JsonResponse
from django.views.generic.base import View
from django.urls import reverse
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.hashers import make_password
from django.db.models import Q

from operations.models import Banner, Advertise, UserInterest, EmailVerifyCode
from projects.models import Project, ProjectCategory
from users.models import UserProfile
from orders.models import UserPayOrder
from utils.SendVerifyCode import send_email_code
# Create your views here.


class IndexView(View):
    '''测试一下'''
    def get(self, request):
        banner_list = Banner.objects.all()
        advertise_list = Advertise.objects.all().order_by('weight')
        category_weight_list = ProjectCategory.objects.all().order_by('-weight')[:3]
        cwl = [p.name for p in category_weight_list]
        other_project_list = Project.objects.exclude(category__name__in=cwl)
        return render(request, 'index.html', {'banner_list': banner_list,
                                              'advertise_list': advertise_list,
                                              'category_weight_list': category_weight_list,
                                              'other_project_list': other_project_list})


class LoginView(View):
    def get(self, request):
        return render(request, 'login.html')

    def post(self, request):
        username = request.POST.get('username', '')
        password = request.POST.get('password', '')
        if not all([username, password]):
            return JsonResponse({'res': 1})
        user = authenticate(username=username, password=password)
        if not user:
            return JsonResponse({'res': 2})
        else:
            login(request, user)
            return JsonResponse({'res': 0})


class RegisterView(View):
    def get(self, request):
        return render(request, 'reg.html')

    def post(self, request):
        username = request.POST.get('username', '')
        password = request.POST.get('password', '')
        email = request.POST.get('email', '')
        account_type = request.POST.get('account_type', '')
        if not all([username, password, email, account_type]):
            return JsonResponse({'res': 1})
        try:
            UserProfile.objects.create(username=username, password=make_password(password), email=email, account_type=account_type)
            user = UserProfile.objects.get(username=username)
            user.is_active = False
            send_email_code(email, send_type='register', user_id=user.id)
        except Exception as e:
            print(e)
            return JsonResponse({'res': 2})
        else:
            return JsonResponse({'res': 0})


class ActiveView(View):
    def get(self, request, active_code):
        active_code, user_id = active_code.split('LyZcW')
        email_record = EmailVerifyCode.objects.filter(verify_code=active_code)
        if email_record:
            email = email_record[0].email
            user = UserProfile.objects.get(email=email, id=user_id)
            user.is_active = True
            user.save()
            return redirect(reverse('users:login'))
        return redirect(reverse('index'))


class LogoutView(View):
    def get(self, request):
        logout(request)
        return redirect(reverse('index'))


class UserAcctTypeView(View):
    def get(self, request):
        return render(request, 'accttype.html')


class UserCenterView(View):
    def get(self, request):
        return render(request, 'member.html')


class UserCrowdFundingView(View):
    def get(self, request):
        user = request.user
        user_orders = UserPayOrder.objects.filter(user_id=user.id)
        interest_projects = UserInterest.objects.filter(user_id=user.id)
        sponsored_projects = Project.objects.filter(admin_id=user.id)
        return render(request, 'minecrowdfunding.html', {'user_orders': user_orders,
                                                         'interest_projects': interest_projects,
                                                         'sponsored_projects': sponsored_projects, })
