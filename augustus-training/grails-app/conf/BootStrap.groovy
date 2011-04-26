import org.codehaus.groovy.grails.commons.ApplicationHolder

class BootStrap {
    def init = { servletContext ->
        // workaround for GRAILS-4580
        ApplicationHolder.application.domainClasses.each { dc ->
            dc.clazz.count()
        }
    }
     def destroy = {
     }
} 
