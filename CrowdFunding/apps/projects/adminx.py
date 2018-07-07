import xadmin

from projects.models import Project, ProjectItem, Company, ProjectTags, ProjectCategory


xadmin.site.register(Company)
xadmin.site.register(ProjectItem)
xadmin.site.register(Project)
xadmin.site.register(ProjectCategory)
xadmin.site.register(ProjectTags)

